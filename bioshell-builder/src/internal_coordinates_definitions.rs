use std::io;
use std::fs::{self};
use std::collections::HashMap;
use log::{info, error};
use bioshell_cif::{CifData, read_cif_buffer};
use bioshell_io::{open_file, split_into_strings};
use crate::BuilderError;
use crate::BuilderError::InternalAtomDefinitionError;

/// Defines which residue an atom used by ``InternalAtomDefinition`` comes from
///
/// A RelaviveResidueLocator value is used to say whether an atom used to define internal coordinates
/// comes from the very same residue as the one being reconstructed, or maybe from the previous residue.
/// E.g. the definition of an N atom defined by the following line (in the CIF format, see
/// [`InternalCoordinatesDatabase`](InternalCoordinatesDatabase) documentation for an example entry):
/// ``` text
/// "ALA prev ' N  ' prev ' CA ' prev ' C  ' this ' N  ' N  1.328685 114.0  180.0 Psi"
/// ```
/// is based on the `N` atom of the ``previous`` residue as well as on the `CA` and `C` atoms of this residue.
#[derive(Clone)]
pub enum RelaviveResidueLocator {
    /// the atom is located in the residue preceding the reconstructed one
    Previous,
    /// the atom is located in the residue being reconstructed
    This,
    /// the atom is located in the residue following the reconstructed one
    Next
}

impl TryFrom<&str> for RelaviveResidueLocator {
    type Error = BuilderError;

    /// Returns a ``RelaviveResidueLocator`` for its string name.
    ///
    /// Three spelling variants are allowed, as shown in the example below
    /// # Example
    /// ```rust
    /// use bioshell_builder::RelaviveResidueLocator;
    /// assert_eq!(RelaviveResidueLocator::try_from("Next").unwrap(), RelaviveResidueLocator::Next);
    /// assert_eq!(RelaviveResidueLocator::try_from("next").unwrap(), RelaviveResidueLocator::Next);
    /// assert_eq!(RelaviveResidueLocator::try_from("NEXT").unwrap(), RelaviveResidueLocator::Next);
    /// ```
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "This" => Ok(RelaviveResidueLocator::This),
            "Next" => Ok(RelaviveResidueLocator::Next),
            "Prev" => Ok(RelaviveResidueLocator::Previous),
            "this" => Ok(RelaviveResidueLocator::This),
            "next" => Ok(RelaviveResidueLocator::Next),
            "prev" => Ok(RelaviveResidueLocator::Previous),
            "THIS" => Ok(RelaviveResidueLocator::This),
            "NEXT" => Ok(RelaviveResidueLocator::Next),
            "PREV" => Ok(RelaviveResidueLocator::Previous),
            _ => {Err(InternalAtomDefinitionError{ error: "Can't find a RelativeResidueLocator for the string".to_string() })
            }
        }
    }
}

impl TryFrom<&RelaviveResidueLocator> for i8 {
    type Error = BuilderError;

    /// Provides integer offset for an enum value.
    ///
    /// This method returns `-1`, `0` or `1` for the `Previous`, `This` and `Next` residue, respectively
    fn try_from(value: &RelaviveResidueLocator) -> Result<Self, Self::Error> {
        match value {
            RelaviveResidueLocator::This => Ok(0),
            RelaviveResidueLocator::Next => Ok(1),
            RelaviveResidueLocator::Previous => Ok(-1),
        }
    }
}

/// Defines an atom by providing it's internal coordinates and the reference frame.
///
/// The three atoms: `a`, `b` and `c` are used to construct a reference frame for an atom `d`.
/// The position of the `d` atom is defined by three internal coordinates:
///    * `r` - distance between `c` and `d`
///    * `planar` - planar angle between `b`, `c` and `d`
///    * `dihedral` - dihedral angle between `a`, `b`, `c` and `d`
/// ```text
///         c----d
///        /
///       /
/// a----b
/// ```
/// The three atoms: `a`, `b` and `c` are referred by their names and ``RelaviveResidueLocator`` enums.
/// Dihedral angle should be defined according to the IUPAC definition, which states that:
/// > *In a Newman projection the torsion angle is the angle (having an absolute value between 0° and 180°)
/// > between bonds to two specified (fiducial) groups, one from the atom nearer (proximal) to the observer
/// > and the other from the further (distal) atom. The torsion angle between groups A and D is then considered
/// > to be positive if the bond A-B is rotated in a clockwise direction through less than 180° in order that
/// > it may eclipse the bond C-D: a negative torsion angle requires rotation in the opposite sense.*
///
/// Each of the four [`RelaviveResidueLocator`](RelaviveResidueLocator) parameters can take values
/// [`This`](RelaviveResidueLocator::This), [`Previous`](RelaviveResidueLocator::Previous) or [`Next`](RelaviveResidueLocator::Next),
/// when an atom belongs to the residue being reconstructed, the preceding or the following one, respectively.
///
#[derive(Clone)]
pub struct InternalAtomDefinition {
    /// Name of a residue this atom belongs to
    pub res_name: String,
    /// Name of the atom this struct defines, i.e. the name of the atom `d`
    pub name: String,
    /// Name of the atomic element for the atom `d`
    pub element: String,
    /// Name of the atom at `a` position
    pub a_name: String,
    /// Name of the atom at `b` position
    pub b_name: String,
    /// Name of the atom at `c` position
    pub c_name: String,
    /// relative index defining the residue the `a` atom belongs to
    pub a_residue: RelaviveResidueLocator,
    /// relative index defining the residue the `b` atom belongs to
    pub b_residue: RelaviveResidueLocator,
    /// relative index defining the residue the `c` atom belongs to
    pub c_residue: RelaviveResidueLocator,
    /// relative index defining the residue the `d` atom belongs to
    pub d_residue: RelaviveResidueLocator,
    /// the distance between `c` and `d` atoms
    pub r: f64,
    /// the planar angle between `b`, `c` and `d` atoms
    pub planar: f64,
    /// the dihedral angle between `a`, `b`, `c` and `d` atoms
    pub dihedral: f64,
    /// Name of this dihedral name, such as `Phi`, `Psi` or `Chi3`
    pub dihedral_name: String
}

impl InternalAtomDefinition {
    /// Creates a new  [`InternalAtomDefinition`](InternalAtomDefinition) struct by filling all its fields
    pub fn from_properties(res_name: &str, name: &str, element: &str, a_locator: RelaviveResidueLocator, a_name: &str,
            b_locator: RelaviveResidueLocator, b_name: &str,
            c_locator: RelaviveResidueLocator, c_name: &str, d_locator: RelaviveResidueLocator,
            r: f64, planar_radians: f64, dihedral_radians: f64, dihedral_name: &str) -> InternalAtomDefinition {

        return InternalAtomDefinition{
            res_name: res_name.to_string(),
            name: name.to_string(), element: element.to_string(),
            a_name: a_name.to_string(),
            b_name: b_name.to_string(), c_name: c_name.to_string(),
            a_residue: a_locator,
            b_residue: b_locator,
            c_residue: c_locator,
            d_residue: d_locator,
            r, planar: planar_radians, dihedral: dihedral_radians,
            dihedral_name: dihedral_name.to_string()
        };
    }

    /// Creates an [`InternalAtomDefinition`](InternalAtomDefinition) struct from a single line in the CIF format.
    ///
    /// Note: data entries must be provided in the proper order!
    /// # Examples
    ///
    /// ```rust
    /// use bioshell_builder::InternalAtomDefinition;
    /// let def = InternalAtomDefinition::from_cif_line("'ALA' this ' N  ' this ' CA ' this ' C  ' next ' N  '  N  1.328685 114.0  180.0 psi");
    /// let def_ok = def.ok().unwrap();
    /// assert_eq!(def_ok.a_name, " N  ".to_string());
    /// assert_eq!(def_ok.name, " N  ".to_string());
    /// ```
    pub fn from_cif_line(line: &str) -> Result<InternalAtomDefinition, BuilderError> {
        let tokens: Vec<String> = split_into_strings(line, true);
        return InternalAtomDefinition::from_strings(&tokens);
    }

    /// Creates an [`InternalAtomDefinition`](InternalAtomDefinition) struct from its properties.
    ///
    /// This method assumes all values to be Strings and parses them accordingly
    ///
    /// # Examples
    ///
    /// ```rust
    /// use bioshell_builder::InternalAtomDefinition;
    /// use bioshell_io::split_into_strings;
    /// let tokens = split_into_strings("'ALA' this ' N  ' this ' CA ' this ' C  ' next ' N  ' N  1.328685 114.0  180.0 psi", false);
    /// let def = InternalAtomDefinition::from_strings(&tokens);
    /// let def_ok = def.ok().unwrap();
    /// assert_eq!(def_ok.a_name, " N  ".to_string());
    /// assert_eq!(def_ok.name, " N  ".to_string());
    /// ```
    pub fn from_strings(tokens: &Vec<String>) -> Result<InternalAtomDefinition, BuilderError> {
        let res_name = &tokens[0];
        let a_residue = RelaviveResidueLocator::try_from(tokens[1].as_str());
        let a_name = &tokens[2];
        let b_residue = RelaviveResidueLocator::try_from(tokens[3].as_str());
        let b_name = &tokens[4];
        let c_residue = RelaviveResidueLocator::try_from(tokens[5].as_str());
        let c_name = &tokens[6];
        let d_residue = RelaviveResidueLocator::try_from(tokens[7].as_str());
        let d_name = &tokens[8];
        let d_element = &tokens[9];
        let r = match tokens[10].parse::<f64>() {
            Ok(r) => r,
            Err(_) => return Err(InternalAtomDefinitionError { error: "Can't parse bond length".to_string() })
        };
        let planar = match tokens[11].parse::<f64>() {
            Ok(p) => p.to_radians(),
            Err(_) => return Err(InternalAtomDefinitionError { error: "Can't parse planar angle".to_string() })
        };
        let dihedral = match tokens[12].parse::<f64>() {
            Ok(d) => d.to_radians(),
            Err(_) => return Err(InternalAtomDefinitionError { error: "Can't parse dihedral angle".to_string() })
        };
        let dihedral_name = &tokens[13];
        return Ok(InternalAtomDefinition{
            res_name: unquoted_substr(res_name).to_string(),
            name: unquoted_substr(d_name).to_string(),
            element: d_element.to_string(),
            a_name: unquoted_substr(a_name).to_string(),
            b_name: unquoted_substr(b_name).to_string(),
            c_name: unquoted_substr(c_name).to_string(),
            a_residue: a_residue?,
            b_residue: b_residue?,
            c_residue: c_residue?,
            d_residue: d_residue?,
            r, planar, dihedral,
            dihedral_name: dihedral_name.clone()
        });
    }
}

/// Removes quotation marks from either end of a string.
///
/// Any other possible quotation marks are left intact
fn unquoted_substr(s:&str) -> &str {
    let from = if s.starts_with('\'') || s.starts_with('"') {1} else {0} as usize;
    let to = if s.ends_with('\'') || s.ends_with('"') {s.len()-1} else {s.len()};

    return &s[from..to];
}

/// Knows how to build a residue or a molecule from it's internal definition.
///
/// # Example
/// ```rust
/// use std::io::BufReader;
/// use bioshell_builder::InternalCoordinatesDatabase;
/// use bioshell_cif::read_cif_buffer;
/// // --- The following CIF-formatted text provides an example definition format
/// const GLY_HEAVY_CIF: &str = "data_GLY
/// loop_
/// _res_name
/// _atom_a_residue_locator
/// _atom_a_name
/// _atom_b_residue_locator
/// _atom_b_name
/// _atom_c_residue_locator
/// _atom_c_name
/// _atom_d_residue_locator
/// _atom_d_name
/// _atom_d_element
/// _c_d_bond_length
/// _b_c_d_planar_angle
/// _a_b_c_d_dihedral_angle
/// _dihedral_angle_name
/// 'GLY' prev ' N  ' prev ' CA ' prev ' C  ' this ' N  ' N  1.328685 114.0  180.0 psi
/// 'GLY' prev ' CA ' prev ' C  ' this ' N  ' this ' CA ' C  1.458001 123.0  180.0 omega
/// 'GLY' prev ' C  ' this ' N  ' this ' CA ' this ' C  ' C  1.523258 110.0 -180.0 phi
/// 'GLY' next ' N  ' this ' CA ' this ' C  ' this ' O  ' O  1.231015 121.0  180.0 -
/// #";
/// // --- Read a CIF block and parse it
/// let mut cif_reader = BufReader::new(GLY_HEAVY_CIF.as_bytes());
/// let data_blocks = read_cif_buffer(&mut cif_reader);
/// // --- Initialize empty database and upload the GLY monomer definition from a CIF block
/// let mut idb = InternalCoordinatesDatabase::new();
/// idb.load_from_cif_data(data_blocks);
/// // --- Now it's possible to fetch the definition of GLY residue from a database
/// let gly_def = idb.get_definition("GLY");
/// assert!(gly_def.is_some());
/// let gly_def = gly_def.unwrap();
/// assert_eq!(gly_def.len(), 4);
/// ```
pub struct InternalCoordinatesDatabase {
    map: HashMap<String, Vec<InternalAtomDefinition>>
}

impl InternalCoordinatesDatabase {

    /// Creates a new, empty database
    pub fn new() -> InternalCoordinatesDatabase { InternalCoordinatesDatabase{ map: Default::default() } }

    /// Creates a database and loads all entries defined by CIF files found in a given folder
    pub fn from_cif_directory(folder: &str) -> Result<InternalCoordinatesDatabase, io::Error> {

        let mut out = InternalCoordinatesDatabase::new();
        info!("Looking for CIF monomer files in: {:?}", &folder);
        for entry in fs::read_dir(folder)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_file() {
                if let Some(ext) =  path.extension() {
                    if ext == "cif" || ext == "CIF" {
                        let fname = path.to_str().unwrap();
                        info!("Reading a file in CIF format: {}", &fname);
                        let mut cif_reader = open_file(fname)?;
                        let data_blocks = read_cif_buffer(&mut cif_reader);
                        out.load_from_cif_data(data_blocks);
                    }
                }
            }
        }

        return Ok(out);
    }

    /// Loads all residue definitions from a buffer providing CIF-formatted data
    pub fn load_from_cif_data(&mut self, data: Vec<CifData>) {
        for block in data {
            if let Some(loop_block) = block.loop_blocks().next() {
                let name = block.name();
                let vec = self.map.entry(name.to_string()).or_insert(vec![]);
                for row in loop_block.rows() {
                    let def = InternalAtomDefinition::from_strings(row);
                    match def {
                        Ok(block) => { vec.push(block) }
                        Err(err) => { error!("{}",err.to_string()) }
                    };
                }
            }
        }
    }

    /// Provides a definition for all atoms of a requested monomer
    pub fn get_definition(&self, residue_name: &str) -> Option<&Vec<InternalAtomDefinition>> {
        return self.map.get(residue_name);
    }

    /// Counts the residue type definitions know to this database
    pub fn count_definitions(&self) -> usize { self.map.len() }
}