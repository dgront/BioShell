use std::cmp::Ordering;
use std::fmt::{Display, Formatter};
use std::string::String;
use crate::calc::Vec3;
use crate::SecondaryStructureTypes;

/// Atom record as found in a single line of a PDB file.
///
///The struct holds all data parsed from [`ATOM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)
/// or [`HETATM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM) lines.
///
/// # Examples
///```rust
/// use bioshell_pdb::PdbAtom;
/// let pdb_line = "ATOM    320  CA  PHE A  43      16.101   9.057  19.587  1.00 18.18           C  ";
/// let a = PdbAtom::from_atom_line(pdb_line);
/// assert_eq!(a.name.as_str(), " CA ");
/// assert_eq!(a.res_name.as_str(), "PHE");
/// assert_eq!(a.is_hetero_atom, false);
/// assert_eq!(format!("{}",a), pdb_line);
///```
#[derive(Clone, Debug)]
pub struct PdbAtom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: String,
    pub entity_id: String,
    pub res_seq: i32,
    pub i_code: char,
    pub pos: Vec3,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: Option<String>,
    pub charge: Option<String>,
    pub is_hetero_atom: bool,
    pub secondary_struct_type: u8,
}

impl PdbAtom {
    /// Returns a default atom.
    ///
    /// By default, an atom is set to alpha-carbon of `ALA1` residue in chain "A", located at `[0,0,0]`
    pub fn new() -> PdbAtom {
        PdbAtom {
            serial: 1,
            name: String::from(" CA "),
            alt_loc: ' ',
            res_name: String::from("ALA"),
            chain_id: String::from("A"),
            entity_id: String::from("1"),
            res_seq: 1,
            i_code: ' ',
            pos: Vec3::from_float(0.0),
            occupancy: 1.0,
            temp_factor: 0.0,
            element: Some(String::from("C")),
            is_hetero_atom: false,
            secondary_struct_type: SecondaryStructureTypes::Coil as u8,
            charge: None
        }
    }

    /// Creates a [`PdbAtom`] by parsing an `ATOM` or `HETATM` record.
    ///
    /// The method automatically sets the [`PdbAtom::is_hetero_atom`](PdbAtom::is_hetero_atom) flag
    /// based on the record type.
    /// The format of the records is defined on the [www.wwpdb.org](https://www.wwpdb.org/) site:
    /// [`ATOM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)
    /// and [`HETATM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM), respectively
    ///
    /// ```
    /// use bioshell_pdb::PdbAtom;
    /// let a = PdbAtom::from_atom_line("ATOM     33  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C");
    /// assert_eq!(a.res_seq, -3);
    /// assert_eq!(a.alt_loc, 'A');
    /// assert_eq!(a.name, " CA ");
    /// assert_eq!(a.element.unwrap(), "C");
    /// let a = PdbAtom::from_atom_line("ATOM     33  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67");
    /// assert_eq!(a.element, None);
    /// ```
    pub fn from_atom_line(pdb_line: &str) -> PdbAtom {
        let serial = pdb_line[6..11].trim().parse::<i32>().unwrap();
        let name = pdb_line[12..16].to_string();
        let alt_loc = pdb_line[16..17].chars().next().unwrap();
        let res_name = pdb_line[17..20].to_string();
        let chain_id = pdb_line[21..22].to_string();
        let res_seq = pdb_line[22..26].trim().parse::<i32>().unwrap();
        let i_code = pdb_line[26..27].chars().next().unwrap();
        let pos = Vec3::from_pdb_line(pdb_line);
        let occupancy = pdb_line[54..60].trim().parse::<f64>().unwrap();
        let temp_factor = pdb_line[60..66].trim().parse::<f64>().unwrap();
        let element = if pdb_line.len()>=78 { Some(pdb_line[77..].trim().to_string()) } else { None };
        return PdbAtom {
            serial,
            name,
            alt_loc,
            res_name,
            chain_id: chain_id.clone(),
            entity_id: chain_id,
            res_seq,
            i_code,
            pos,
            occupancy,
            temp_factor,
            element,
            charge: None,
            is_hetero_atom: pdb_line.starts_with("H"),
            secondary_struct_type: 12
        };
    }

    pub fn hec_code(&self) -> u8 {
        let sec_str_type = SecondaryStructureTypes::from_pdb_class(self.secondary_struct_type as usize);
        SecondaryStructureTypes::hec_code(&sec_str_type)
    }
}

impl PartialEq<Self> for PdbAtom {
    /// Returns true if two [`PdbAtom`](PdbAtom)s are equal.
    ///
    /// The equality of [`PdbAtom`](PdbAtom)s implies that:
    ///   - they belong to the same chain
    ///   - they belong to the same residue, as identified by their `res_seq` and `i_code` fields
    ///   - their serial numbers are identical
    fn eq(&self, other: &Self) -> bool {

        self.chain_id == other.chain_id && self.res_seq==other.res_seq
            && self.i_code==self.i_code && self.serial==self.serial
    }
}

impl PartialOrd<Self> for PdbAtom {

    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.chain_id < other.chain_id { return Some(Ordering::Less) }
        if self.chain_id > other.chain_id { return Some(Ordering::Greater) }
        if self.res_seq < other.res_seq { return Some(Ordering::Less) }
        if self.res_seq > other.res_seq { return Some(Ordering::Greater) }
        if self.i_code < other.i_code { return Some(Ordering::Less) }
        if self.i_code > other.i_code { return Some(Ordering::Greater) }
        if self.serial < other.serial { return Some(Ordering::Less) }
        if self.serial > other.serial { return Some(Ordering::Greater) }

        return  Some(Ordering::Equal);
    }
}

impl Eq for PdbAtom {}

impl Ord for PdbAtom {
    fn cmp(&self, other: &Self) -> Ordering { return self.partial_cmp(&other).unwrap(); }
}

impl Display for PdbAtom {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let chrg = match &self.charge {
            None => {"  "}
            Some(c) => {c}
        };
        let elem = match &self.element {
            None => {"  "}
            Some(c) => {c}
        };
        if self.is_hetero_atom {
            write!(f, "HETATM{:>5} {:^4}{:>1}{:>3} {:>1}{:>4}{:1}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}{:2}",
                   self.serial, self.name, self.alt_loc, self.res_name, self.chain_id, self.res_seq,
                   self.i_code, self.pos.x, self.pos.y, self.pos.z, self.occupancy, self.temp_factor,
                   elem, chrg)
        } else {
            write!(f, "ATOM  {:>5} {:^4}{:>1}{:>3} {:>1}{:>4}{:1}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}{:2}",
                   self.serial, self.name, self.alt_loc, self.res_name, self.chain_id, self.res_seq,
                   self.i_code, self.pos.x, self.pos.y, self.pos.z, self.occupancy, self.temp_factor,
                   elem, chrg)
        }
    }
}
/// Returns `true` if two given atoms belong to the very same residue
///
/// # Examples
/// ```
/// use bioshell_pdb::{PdbAtom, same_residue_atoms};
/// let a1 = PdbAtom::from_atom_line("ATOM    389  CG2 VAL A  50       7.150   8.278  10.760  1.00 20.57           C");
/// let a2 = PdbAtom::from_atom_line("ATOM    390  N   LEU A  51      10.919   9.836  12.777  1.00 10.30           N");
/// let a3 = PdbAtom::from_atom_line("ATOM    391  CA  LEU A  51      12.088   9.803  13.653  1.00  9.53           C  ");
/// assert!(!same_residue_atoms(&a1, &a2));
/// assert!(same_residue_atoms(&a2, &a3));
/// ```
pub fn same_residue_atoms(ai: &PdbAtom, aj: &PdbAtom) -> bool {
    ai.res_seq==aj.res_seq && ai.i_code==aj.i_code
}