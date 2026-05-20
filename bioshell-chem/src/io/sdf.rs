use std::io::{BufRead, Lines};
use crate::{Atom, BondType, ChemErrors, Element, Molecule};
use crate::io::{parse, parse_substr};

/// Reads a molecule from an SDF file.
///
/// The SDF format is described [on Wikipedia](https://en.wikipedia.org/wiki/Chemical_table_file#SDF).
/// A simple example (ethanol molecule) is given below:
/// ```text
#[doc = include_str!("../../tests/test_files/EOH_small.sdf")]
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
        let a = parse_atom_line(atom_idx, &line)?;
        mol.add_atom(a)?;
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
    let n_atoms: usize = parse_substr(line, 0, 3, "SDF atom count")?;
    let n_bonds: usize = parse_substr(line, 3, 6, "SDF bond count")?;

    Ok((n_atoms, n_bonds))
}

fn parse_atom_line(idx: usize, line: &str) -> Result<Atom, ChemErrors> {
    let fields: Vec<&str> = line.split_whitespace().collect();

    if fields.len() < 4 {
        return Err(ChemErrors::IncorrectSdfFormat(format!("Invalid SDF atom line: {line}")));
    }

    let vx: f64 = parse(fields[0], "sdf atom x")?;
    let vy: f64 = parse(fields[1], "sdf atom y")?;
    let vz: f64 = parse(fields[2], "sdf atom z")?;

    let e = parse_atom_element(line)?;
    let q = parse_atom_charge(line)?;

    let mut a = Atom::charged(idx, e, q);
    a.set_pos3(vx, vy, vz);

    Ok(a)
}

fn parse_atom_charge(line: &str) -> Result<i8, ChemErrors> {
    let charge_code: i8 =parse_substr(line, 36, 39, "SDF atom charge field")?;

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
        .ok_or_else(|| ChemErrors::LineTooShort(line.into(), 34))?
        .trim()
        .parse()
}

fn parse_bond_line(line: &str) -> Result<(usize, usize, BondType), ChemErrors> {
    let a: usize = parse_substr(line, 0, 3, "SDF bond from")?;
    let b: usize = parse_substr(line, 3, 6, "SDF bond to")?;
    let order: u8 = parse_substr(line, 6, 9, "SDF bond order")?;

    let bond_type = match order {
        1 => BondType::Single,
        2 => BondType::Double,
        3 => BondType::Triple,
        4 => BondType::Aromatic,
        _ => BondType::Unknown
    };

    Ok((a, b, bond_type))
}