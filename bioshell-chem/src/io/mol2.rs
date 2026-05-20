use std::io::{BufRead, BufReader, Read};
use crate::{Atom, ChemErrors, Molecule, BondType, Element};
use crate::ChemErrors::{InvalidMol2AtomId, InvalidMol2AtomLine, InvalidMol2BondLine, NumericParsingError, UnknownElement};


/// Reads a molecule from a MOL2 file.
///
/// A file in `.mol2` format lists atoms and bonds of  a single molecule. An example is given below:
/// ```text
#[doc = include_str!("../../tests/test_files/EOH.mol2")]
///```
///
/// Such a file can be loaded as:
/// ```
/// use bioshell_chem::{ChemErrors};
/// # fn main() -> Result<(), ChemErrors> {
/// use bioshell_core::io::open_file;
/// use bioshell_chem::molecule_from_sdf;
/// let reader = open_file("./tests/test_files/toluene.sdf")?;
/// let mol = molecule_from_sdf(reader)?;
/// assert_eq!(mol.count_atoms(), 15);
/// assert_eq!(mol.count_bonds(), 15);
/// # Ok(())
/// # }
/// ```
pub fn molecule_from_mol2<R: Read>(reader: R) -> Result<Molecule, ChemErrors> {
    let reader = BufReader::new(reader);

    let mut molecule = Molecule::new("");
    let mut section: Option<&str> = None;

    for line_result in reader.lines() {
        let line = line_result?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        if line.starts_with("@<TRIPOS>") {
            section = match line {
                "@<TRIPOS>MOLECULE" => Some("MOLECULE"),
                "@<TRIPOS>ATOM" => Some("ATOM"),
                "@<TRIPOS>BOND" => Some("BOND"),
                _ => None,
            };
            continue;
        }

        match section {
            Some("ATOM") => {
                let fields: Vec<&str> = line.split_whitespace().collect();

                if fields.len() < 6 {
                    return Err(InvalidMol2AtomLine(line.to_string()));
                }

                let atom_id: usize = fields[0]
                    .parse()
                    .map_err(|_| InvalidMol2AtomId(fields[0].to_string()))?;

                let mol2_atom_type = fields[5];
                let element = element_from_mol2_type(mol2_atom_type)?;
                let mut a = Atom::neutral(atom_id - 1, element);
                // let atom_name = fields[1];
                let x: f64 = fields[2].parse().map_err(|_| NumericParsingError("x".into(), fields[2].to_string()))?;
                let y: f64 = fields[3].parse().map_err(|_| NumericParsingError("y".into(), fields[3].to_string()))?;
                let z: f64 = fields[4].parse().map_err(|_| NumericParsingError("z".into(), fields[4].to_string()))?;
                a.set_pos3(x, y, z);
                molecule.add_atom(a)?;
            }

            Some("BOND") => {
                let fields: Vec<&str> = line.split_whitespace().collect();

                if fields.len() < 4 { return Err(InvalidMol2BondLine(line.to_string())); }

                let i: usize = fields[1].parse().map_err(|_| InvalidMol2AtomId(fields[1].to_string()))?;
                let j: usize = fields[2].parse().map_err(|_| InvalidMol2AtomId(fields[2].to_string()))?;
                let bond_type = fields[3];

                molecule.bind_atoms(i - 1, j - 1, BondType::from_code(bond_type))?;
            }

            _ => {}
        }
    }

    Ok(molecule)
}

fn element_from_mol2_type(atom_type: &str) -> Result<Element, ChemErrors> {
    let symbol = atom_type
        .split('.')
        .next()
        .ok_or_else(|| UnknownElement(atom_type.to_string()))?;

    symbol.parse().map_err(|_| UnknownElement(atom_type.to_string()))
}