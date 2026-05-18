use std::io::{BufRead, Lines};
use crate::{Atom, BondType, ChemErrors, Element, Molecule};


/// Reads a molecule from an SDF file.
///
/// The SDF format is described here: <https://en.wikipedia.org/wiki/Chemical_table_file#SDF>
///
/// # Example
/// ```
/// use bioshell_chem::{ChemErrors};
/// # fn main() -> Result<(), ChemErrors> {
/// use bioshell_core::io::open_file;
/// use bioshell_chem::molecule_from_sdf;
/// let reader = open_file("./tests/test_files/dspc.sdf")?;
/// let mol = molecule_from_sdf(reader)?;
/// assert_eq!(mol.count_atoms(), 142);
/// assert_eq!(mol.count_bonds(), 141);
/// # Ok(())
/// # }
/// ```
pub fn molecule_from_sdf<R: BufRead>(reader: R) -> Result<Molecule, ChemErrors> {
    let mut lines = reader.lines();

    let name = lines.next_sdf_line("Missing SDF title line")?;

    lines.next_sdf_line("Missing SDF program line")?;
    lines.next_sdf_line("Missing SDF comment line")?;

    let counts = lines.next_sdf_line("Missing SDF counts line")?;
    let (n_atoms, n_bonds) = parse_counts_line(&counts)?;

    let mut mol = Molecule::new(&name);

    for atom_idx in 0..n_atoms {
        let line = lines.next_sdf_line("Missing SDF atom line")?;
        let element = parse_atom_element(&line)?;
        let charge = parse_atom_charge(&line)?;
        mol.add_atom(Atom::charged(atom_idx, element, charge))?;
    }

    for _ in 0..n_bonds {
        let line = lines.next_sdf_line("Missing SDF bond line")?;
        let (a, b, bond_type) = parse_bond_line(&line)?;
        mol.bind_atoms(a - 1, b - 1, bond_type)?;
    }

    Ok(mol)
}

// ---------- Helper functions for reading SFD files

trait SdfLineExt {
    fn next_sdf_line(&mut self, msg: &'static str) -> Result<String, ChemErrors>;
}

impl<R: BufRead> SdfLineExt for Lines<R> {
    fn next_sdf_line(&mut self, msg: &'static str) -> Result<String, ChemErrors> {
        self.next()
            .transpose()?
            .ok_or_else(|| ChemErrors::IncorrectSdfFormat(msg.into()))
    }
}

/// Extracts the number of atoms and the number of bonds from a line
fn parse_counts_line(line: &str) -> Result<(usize, usize), ChemErrors> {
    let n_atoms = line.get(0..3)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF counts line".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid atom count".into()))?;

    let n_bonds = line.get(3..6)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF counts line".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid bond count".into()))?;

    Ok((n_atoms, n_bonds))
}

fn parse_atom_charge(line: &str) -> Result<i8, ChemErrors> {
    let charge_code: i8 = line.get(36..39)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF atom charge field".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid SDF atom charge code".into()))?;

    match charge_code {
        0 => Ok(0),
        1 => Ok(3),
        2 => Ok(2),
        3 => Ok(1),
        4 => Ok(0), // doublet radical; not formal charge
        5 => Ok(-1),
        6 => Ok(-2),
        7 => Ok(-3),
        _ => Err(ChemErrors::IncorrectSdfFormat(
            format!("Unsupported SDF atom charge code: {charge_code}")
        )),
    }
}

fn parse_atom_element(line: &str) -> Result<Element, ChemErrors> {
    line.get(31..34)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF atom line".into()))?
        .trim()
        .parse()
}

fn parse_bond_line(line: &str) -> Result<(usize, usize, BondType), ChemErrors> {
    let a: usize = line
        .get(0..3)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF bond line".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid first bond atom".into()))?;

    let b: usize = line
        .get(3..6)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF bond line".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid second bond atom".into()))?;

    let order: u8 = line
        .get(6..9)
        .ok_or_else(|| ChemErrors::IncorrectSdfFormat("Invalid SDF bond line".into()))?
        .trim()
        .parse()
        .map_err(|_| ChemErrors::IncorrectSdfFormat("Invalid SDF bond order".into()))?;

    let bond_type = match order {
        1 => BondType::Single,
        2 => BondType::Double,
        3 => BondType::Triple,
        4 => BondType::Aromatic,
        _ => {
            return Err(ChemErrors::IncorrectSdfFormat(format!("Unsupported SDF bond order: {order}")))
        }
    };

    Ok((a, b, bond_type))
}