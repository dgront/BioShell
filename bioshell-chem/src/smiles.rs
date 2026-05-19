use crate::{Atom, BondType, ChemErrors, Element, Molecule};
use crate::ChemErrors::IncorrectSmilesString;



impl Molecule {
    pub fn from_smiles(name: &str, s: &str) -> Result<Self, ChemErrors> {
        let tokens = Self::tokenize(s)?;
        Self::from_smiles_tokens(name, &tokens)
    }

    fn tokenize(smiles: &str) -> Result<Vec<&str>, ChemErrors> {
        let mut tokens = Vec::new();
        let mut i = 0;

        while i < smiles.len() {
            let rest = &smiles[i..];
            let ch = rest.chars().next().unwrap();

            let token_len = match ch {
                '[' => {
                    let end = rest.find(']').ok_or_else(|| {
                        IncorrectSmilesString{ smiles: smiles.to_string(), reason: "Unclosed bracket atom in SMILES".into() }
                    })?;
                    end + 1
                }

                '%' => {
                    if rest.len() < 3 || !rest[1..3].chars().all(|c| c.is_ascii_digit()) {
                        return Err(IncorrectSmilesString{ smiles: smiles.to_string(), reason: "Invalid SMILES ring index after '%'".into() });
                    }
                    3
                }

                'A'..='Z' => {
                    let mut len = ch.len_utf8();

                    if let Some(next) = rest[len..].chars().next() {
                        if next.is_ascii_lowercase() {
                            len += next.len_utf8();
                        }
                    }

                    len
                }

                'a'..='z' => ch.len_utf8(),

                '0'..='9' => ch.len_utf8(),

                '-' | '=' | '#' | ':' | '(' | ')' | '.' | '/' | '\\' => ch.len_utf8(),

                _ => {
                    return Err(IncorrectSmilesString{ smiles: smiles.to_string(), reason: format!("Unsupported SMILES character: {ch}").to_string() });
                }
            };

            tokens.push(&smiles[i..i + token_len]);
            i += token_len;
        }

        Ok(tokens)
    }


    pub fn from_smiles_tokens(name: &str, tokens: &[&str]) -> Result<Molecule, ChemErrors> {
        let mut mol = Molecule::new(name);

        let mut current_atom: Option<usize> = None;
        let mut pending_bond = BondType::Single;
        let mut branch_stack: Vec<usize> = Vec::new();
        let mut ring_openings: std::collections::HashMap<String, (usize, BondType)> =
            std::collections::HashMap::new();

        for token in tokens {
            match *token {
                // --- Bonds ---
                "-" => pending_bond = BondType::Single,
                "=" => pending_bond = BondType::Double,
                "#" => pending_bond = BondType::Triple,
                ":" => pending_bond = BondType::Aromatic,

                // --- We started a new branch: push the current atom index to the stack ---
                "(" => {
                    let atom = current_atom.ok_or_else(|| IncorrectSmilesString {
                        smiles: tokens.join(""),
                        reason: "Branch opened before any atom".into(),
                    })?;
                    branch_stack.push(atom);
                }

                // --- Closing the current branch: pop the atom index from the stack ---
                ")" => {
                    current_atom = Some(branch_stack.pop().ok_or_else(|| {
                        IncorrectSmilesString {
                            smiles: tokens.join(""),
                            reason: "Branch closed without matching opening".into(),
                        }
                    })?);
                }
                // --- Dot: end of current molecule, start a new disconnected one ---
                "." => {
                    current_atom = None;
                    pending_bond = BondType::Single; // USE NONE HERE?
                }
                // --- Directional bond markers (for stereochemistry) ---
                "/" | "\\" => {
                    // For now: ignore directional bond markers.
                    // Later they should be stored as stereochemical bond annotations.
                }

                token if Self::is_ring_token(token) => {
                    let atom = current_atom.ok_or_else(|| IncorrectSmilesString {
                        smiles: tokens.join(""),
                        reason: "Ring closure without current atom".into(),
                    })?;

                    if let Some((other_atom, opening_bond)) = ring_openings.remove(token) {
                        let bond = if pending_bond == BondType::Single {
                            opening_bond
                        } else {
                            pending_bond
                        };

                        mol.bind_atoms(other_atom, atom, bond)?;
                    } else {
                        ring_openings.insert(token.to_string(), (atom, pending_bond));
                    }

                    pending_bond = BondType::Single;
                }

                token if Self::is_atom_token(token) => {
                    let (element, charge, aromatic) = Self::parse_smiles_atom_token(token)?;
                    let atom_idx = mol.count_atoms();

                    mol.add_atom(Atom::charged(atom_idx, element, charge))?;
                    if let Some(prev) = current_atom {
                        let bond = if pending_bond == BondType::Single && aromatic {
                            BondType::Aromatic
                        } else {
                            pending_bond
                        };

                        mol.bind_atoms(prev, atom_idx, bond)?;
                    }
                    current_atom = Some(atom_idx);
                    if aromatic {
                        pending_bond = BondType::Aromatic;
                    } else { pending_bond = BondType::Single; }
                }

                _ => {
                    return Err(IncorrectSmilesString {
                        smiles: tokens.join(""),
                        reason: format!("Unsupported SMILES token: {token}"),
                    });
                }
            }
        }

        if !branch_stack.is_empty() {
            return Err(IncorrectSmilesString { smiles: tokens.join(""),
                reason: "Unclosed branch in SMILES".into(),
            });
        }

        if !ring_openings.is_empty() {
            return Err(IncorrectSmilesString { smiles: tokens.join(""),
                reason: "Unclosed ring in SMILES".into(), });
        }

        Ok(mol)
    }

    fn is_atom_token(token: &str) -> bool {
        token.starts_with('[')
            || token.chars().next().is_some_and(|c| c.is_ascii_alphabetic())
    }

    fn is_ring_token(token: &str) -> bool {
        token.chars().all(|c| c.is_ascii_digit())
            || token.starts_with('%')
    }

    fn is_aromatic_token_symbol(token: &str) -> bool {
        token.chars().next().is_some_and(|c| c.is_ascii_lowercase())
    }

    fn parse_smiles_atom_token(token: &str) -> Result<(Element, i8, bool), ChemErrors> {
        if token.starts_with('[') {
            return Self::parse_bracket_atom(token);
        }

        let aromatic = token.chars().next().unwrap().is_ascii_lowercase();

        let symbol = if aromatic {
            let mut chars = token.chars();
            let first = chars.next().unwrap().to_ascii_uppercase();
            first.to_string()
        } else {
            token.to_string()
        };

        let element: Element = symbol.parse()?;
        Ok((element, 0, aromatic))
    }

    fn parse_bracket_atom(token: &str) -> Result<(Element, i8, bool), ChemErrors> {
        let inner = token.strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .ok_or_else(|| IncorrectSmilesString {
                smiles: token.to_string(),
                reason: "Invalid bracket atom".into(),
            })?;

        let aromatic = inner.chars().next().is_some_and(|c| c.is_ascii_lowercase());

        let symbol_len = if inner.len() >= 2 {
            let mut chars = inner.chars();
            let first = chars.next().unwrap();
            let second = chars.next().unwrap();

            if first.is_ascii_uppercase() && second.is_ascii_lowercase() {
                2
            } else {
                1
            }
        } else {
            1
        };

        let raw_symbol = &inner[..symbol_len];

        let symbol = if aromatic {
            raw_symbol.to_ascii_uppercase()
        } else {
            raw_symbol.to_string()
        };

        let charge = if inner.contains('-') {
            -1
        } else if inner.contains('+') {
            1
        } else {
            0
        };

        let element: Element = symbol.parse()?;
        Ok((element, charge, aromatic))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_tokenizer() -> Result<(), ChemErrors> {
        let cases = [
            ("CCO", vec!["C", "C", "O"]),
            ("CC(Cl)Br", vec!["C", "C", "(", "Cl", ")", "Br"]),
            ("CC(=O)O", vec!["C", "C", "(", "=", "O", ")", "O"]),
            ("c1ccccc1", vec!["c", "1", "c", "c", "c", "c", "c", "1"]),
            ("C1CCCCC1", vec!["C", "1", "C", "C", "C", "C", "C", "1"]),
            ("C%12CCCCC%12", vec!["C", "%12", "C", "C", "C", "C", "C", "%12"]),
            ("[NH4+]", vec!["[NH4+]"]),
            ("[NH4+]CC(=O)[O-]", vec!["[NH4+]", "C", "C", "(", "=", "O", ")", "[O-]"]),
            ("F/C=C\\Cl", vec!["F", "/", "C", "=", "C", "\\", "Cl"]),
            ("CC.CC", vec!["C", "C", ".", "C", "C"]),
            ("c1ncc[nH]1", vec!["c", "1", "n", "c", "c", "[nH]", "1"]),
            ("c1ccc2cc3ccccc3cc2c1", vec!["c", "1", "c", "c", "c", "2", "c", "c", "3", "c", "c", "c", "c", "c", "3", "c", "c", "2", "c", "1"]),
        ];

        for (smiles, expected) in cases {
            let tokens = Molecule::tokenize(smiles)?;
            assert_eq!(tokens, expected, "failed for SMILES: {smiles}");
        }
        Ok(())
    }
}