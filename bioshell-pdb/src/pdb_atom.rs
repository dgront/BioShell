use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::string::String;
use bioshell_cif::{CifData, CifError, CifTable, entry_has_value, parse_item_or_error, value_or_default};
use bioshell_cif::CifError::ItemParsingError;
use crate::calc::Vec3;
use crate::{HasCartesians, PDBError, SecondaryStructureTypes, Structure};
use crate::PDBError::CifParsingError;

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
    pub secondary_struct_type: SecondaryStructureTypes,
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
            secondary_struct_type: SecondaryStructureTypes::Coil,
            charge: None
        }
    }

    /// Creates a [`PdbAtom`] by parsing an `ATOM` or `HETATM` record of a PBD file.
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
            secondary_struct_type: SecondaryStructureTypes::Coil,
        };
    }

    /// Returns a secondary structure code of the residue this atom belongs to.
    ///
    /// # Examples
    /// ```
    /// use bioshell_pdb::{PdbAtom, SecondaryStructureTypes};
    /// let a = PdbAtom::from_atom_line("ATOM     33  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C");
    /// assert_eq!(a.secondary_struct_type, SecondaryStructureTypes::Coil);
    /// assert_eq!(a.hec_code(), b'C');
    /// ```
    pub fn hec_code(&self) -> u8 { self.secondary_struct_type.hec_code() }

    /// Creates atoms from CIF-formatted data
    pub(crate) fn from_cif_data(cif_data_block: &CifData, structure: &mut Structure) -> Result<(), PDBError> {

        let atoms_tokens = CifTable::new(cif_data_block, "_atom_site.",
     ["id", "label_atom_id",
                     "label_alt_id", "label_comp_id", "auth_asym_id", "label_seq_id", "pdbx_PDB_ins_code",
                     "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "type_symbol",
                     "pdbx_PDB_model_num", "auth_seq_id", "label_entity_id"])?;

        let mut model_ids: HashMap<usize, usize> = HashMap::new();       // links model id to model index in the model_coordinates structure
        for tokens in atoms_tokens.iter() {
            let model_id = tokens[13].parse::<usize>().unwrap();  // model_id of the current atom
            if ! model_ids.contains_key(&model_id) {                    // if we have a new model...
                model_ids.insert(model_id, model_ids.len());
                structure.model_coordinates.push(vec![]);
            }
            let mdl_idx = model_ids[&model_id];
            let x = tokens[7].parse::<f64>().unwrap();
            let y = tokens[8].parse::<f64>().unwrap();
            let z = tokens[9].parse::<f64>().unwrap();
            let pos = Vec3::new(x, y, z);

            if model_ids.len() == 1 {
                structure.model_coordinates[mdl_idx].push(pos.clone());
                match PdbAtom::create_pdb_atom(&tokens, pos) {
                    Ok(a) => { structure.atoms.push(a); }
                    Err(e) => {
                        let data_line = tokens.join(" ");
                        return match e {
                            ItemParsingError { item, type_name, details: _ } => {
                                Err(CifParsingError(ItemParsingError {
                                    item, type_name, details: data_line,
                                }))
                            }
                            _ => { Err(CifParsingError(e)) }
                        }
                    }
                }
            } else {
                structure.model_coordinates[mdl_idx].push(pos);
            }
        }
        Ok(())
    }

    /// A helper function to create an atom based on given string tokens
    fn create_pdb_atom(tokens: &[&str; 16], pos: Vec3) -> Result<PdbAtom, CifError> {

        let serial = parse_item_or_error!(tokens[0], i32);
        let name = format!("{:^4}", tokens[1]);
        let alt_loc = value_or_default(tokens[2], ' ');
        let res_name = tokens[3].to_string();
        let chain_id = tokens[4].to_string();
        let entity_id = tokens[15].to_string();
        let res_seq = if entry_has_value(tokens[14]) {
            parse_item_or_error!(tokens[14], i32)
        } else {
            parse_item_or_error!(tokens[5], i32)
        };

        let i_code = value_or_default(tokens[6], ' ');
        let occupancy = parse_item_or_error!(tokens[10], f64);
        let temp_factor = parse_item_or_error!(tokens[11], f64);
        let element = Some(tokens[12].to_string());
        let charge = None;
        let is_hetero_atom = false;
        let secondary_struct_type = SecondaryStructureTypes::Coil;
        let a = PdbAtom {
            serial, name, alt_loc, res_name,
            chain_id, entity_id, res_seq, i_code, pos, occupancy, temp_factor,
            element, charge, is_hetero_atom, secondary_struct_type,
        };

        return Ok(a);
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

impl HasCartesians for PdbAtom {
    /// Provides Cartesian coordinates of this atom
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{assert_delta, PdbAtom, HasCartesians};
    /// let atom_line = "ATOM   2831  OE1BGLN A 294C    -27.117  12.343  28.479  1.00  9.58           O  ";
    /// let atom = PdbAtom::from_atom_line(atom_line);
    /// assert_delta!(atom.position().x, -27.117, 0.001);
    /// ```
    fn position(&self) -> &Vec3 { &self.pos }
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