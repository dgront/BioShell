use std::fmt::{Display, Formatter};
use std::string::String;
use crate::calc::Vec3;

/// Atom record as found in a single line of a PDB file.
///
///The struct holds all data parsed from [`ATOM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM)
/// or [`HETATM`](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#HETATM) lines.
///
/// # Examples
///```rust
/// use bioshell_pdb::PdbAtom;
/// let pdb_line = "ATOM    320  CA  PHE A  43      16.101   9.057  19.587  1.00 18.18           C ";
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
    pub res_seq: i32,
    pub i_code: char,
    pub pos: Vec3,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: Option<String>,
    pub charge: Option<String>,
    pub is_hetero_atom: bool,
    pub secondary_struct_symbol: char,
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
            res_seq: 1,
            i_code: ' ',
            pos: Vec3::from_float(0.0),
            occupancy: 1.0,
            temp_factor: 0.0,
            element: Some(String::from("C")),
            is_hetero_atom: false,
            secondary_struct_symbol: 'C',
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
    /// assert_eq!(a.alt_loc, "A");
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
        let x = pdb_line[30..38].trim().parse::<f64>().unwrap();
        let y = pdb_line[38..46].trim().parse::<f64>().unwrap();
        let z = pdb_line[46..54].trim().parse::<f64>().unwrap();
        let occupancy = pdb_line[54..60].trim().parse::<f64>().unwrap();
        let temp_factor = pdb_line[60..66].trim().parse::<f64>().unwrap();
        let element = if pdb_line.len()>=78 { Some(pdb_line[77..].trim().to_string()) } else { None };
        return PdbAtom{
            serial,
            name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            pos: Vec3::new(x, y, z),
            occupancy,
            temp_factor,
            element: element,
            charge: None,
            is_hetero_atom: pdb_line.starts_with("H"),
            secondary_struct_symbol: 'C'
        };
    }
}

impl PartialEq for PdbAtom {
    fn eq(&self, other: &Self) -> bool {
        self.serial == other.serial
    }
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