use std::io::BufRead;
use crate::{Atom, BondType, ChemErrors, Element, Molecule};
use crate::io::parse;

pub fn molecule_from_itp<R: BufRead>(reader: R) -> Result<Molecule, ChemErrors> {
    let mut molecule = Molecule::new("");

    let mut section: Option<String> = None;
    let mut expect_molecule_name = false;

    for line_result in reader.lines() {
        let line = line_result?;
        let line = line.trim();

        if line.is_empty() || line.starts_with(';') { continue; }

        if line.starts_with('[') && line.ends_with(']') {
            let name = line
                .trim_start_matches('[')
                .trim_end_matches(']')
                .trim()
                .to_ascii_lowercase();

            expect_molecule_name = name == "moleculetype";
            section = Some(name);
            continue;
        }

        // Remove comments (anything after ';') and trim whitespace
        let data = line.split_once(';').map(|(data, _comment)| data).unwrap_or(line).trim();
        // Skip empty lines after removing comments
        if data.is_empty() { continue; }

        match section.as_deref() {
            Some("moleculetype") if expect_molecule_name => {
                let fields: Vec<&str> = data.split_whitespace().collect();

                if fields.len() < 2 {
                    return Err(ChemErrors::IncorrectItpFormat(
                        format!("Invalid [ moleculetype ] line: {line}"),
                    ));
                }

                molecule.molecule_name = fields[0].to_string();

                expect_molecule_name = false;
            }

            Some("atoms") => {
                let fields: Vec<&str> = data.split_whitespace().collect();

                if fields.len() < 8 {
                    return Err(ChemErrors::IncorrectItpFormat(format!("Invalid [ atoms ] line: {line}")));
                }

                let nr: usize = parse(fields[0], "itp atom nr")?;
                let _atom_type = fields[1];
                let _resnr: usize = parse(fields[2], "itp atom resnr")?;
                let _residu = fields[3];
                let _atom_name = fields[4];
                let _cgnr: usize = parse(fields[5], "itp atom cgnr")?;
                let _charge: f64 = parse(fields[6], "itp atom charge")?;
                let mass: f64 = parse(fields[7], "itp atom mass")?;

                let element = element_by_mass(mass);

                let mut atom = Atom::neutral(nr - 1, element);
                atom.set_pos3(0.0, 0.0, 0.0);

                molecule.add_atom(atom)?;
            }

            Some("bonds") => {
                let fields: Vec<&str> = data.split_whitespace().collect();

                if fields.len() < 3 {
                    return Err(ChemErrors::IncorrectItpFormat(format!("Invalid [ bonds ] line: {line}")));
                }

                let ai: usize = parse(fields[0], "itp bond ai")?;
                let aj: usize = parse(fields[1], "itp bond aj")?;
                let _funct: usize = parse(fields[2], "itp bond funct")?;

                molecule.bind_atoms(ai - 1, aj - 1, BondType::Single)?;
            }

            _ => {}
        }
    }

    Ok(molecule)
}

fn element_by_mass(mass: f64) -> Element {

    let mut best_match = Element::C;
    let mut best_diff = f64::MAX;
    for e in Element::ALL {
        let diff = (e.mass() - mass).abs();
        if diff < best_diff {
            best_diff = diff;
            best_match = *e;
        }
    }
    return best_match;
}