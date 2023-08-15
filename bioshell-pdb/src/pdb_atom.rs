use std::string::String;
use bioshell_seq::chemical::ResidueType;
use crate::pdb_parsing_error::ParseError;

/// Atom data as found in a single line of a PDB file
///
///The struct holds all data parsed from ATOM or HETATM lines
///
/// # Examples
///```rust
/// use bioshell_pdb::PdbAtom;
/// let a = PdbAtom::from_atom_line("ATOM    320  CA  PHE A  43      16.101   9.057  19.587  1.00 18.18           C ");
/// assert_eq!(a.name.as_str(), " CA ");
/// assert_eq!(a.res_name.as_str(), "PHE");
/// assert_eq!(a.is_hetero_atom, false);
///```
#[derive(Clone, Debug)]
pub struct PdbAtom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: String,
    pub res_name: String,
    pub chain_id: String,
    pub res_seq: i32,
    pub i_code: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: Option<String>,
    pub charge: Option<String>,
    pub is_hetero_atom: bool,
    pub secondary_struct_symbol: char,
}

impl PdbAtom {
    pub fn new() -> PdbAtom {
        PdbAtom {
            serial: 1,
            name: String::new(),
            alt_loc: String::new(),
            res_name: String::from("ALA"),
            chain_id: String::new(),
            res_seq: 1,
            i_code: String::new(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 0.0,
            temp_factor: 0.0,
            element: None,
            is_hetero_atom: false,
            secondary_struct_symbol: 'C',
            charge: None
        }
    }

    // pub fn to_string(&self) -> String{
    //     return PdbLineParser::assemble_atom(&self);
    // }

    pub fn from_atom_line(pdb_line: &str) -> PdbAtom {
        let serial = pdb_line[6..11].trim().parse::<i32>().unwrap();
        let name = pdb_line[12..16].to_string();
        let alt_loc = pdb_line[16..17].to_string();
        let res_name = pdb_line[17..20].to_string();
        let chain_id = pdb_line[21..22].to_string();
        let res_seq = pdb_line[22..26].trim().parse::<i32>().unwrap();
        let i_code = pdb_line[26..27].to_string();
        let x = pdb_line[30..38].trim().parse::<f64>().unwrap();
        let y = pdb_line[38..46].trim().parse::<f64>().unwrap();
        let z = pdb_line[46..54].trim().parse::<f64>().unwrap();
        let occupancy = pdb_line[54..60].trim().parse::<f64>().unwrap();
        let temp_factor = pdb_line[60..66].trim().parse::<f64>().unwrap();
        return PdbAtom{
            serial,
            name,
            alt_loc,
            res_name,
            chain_id,
            res_seq,
            i_code,
            x,
            y,
            z,
            occupancy,
            temp_factor,
            element: None,
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