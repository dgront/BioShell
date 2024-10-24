use std::collections::HashMap;
use std::ops::Range;
use crate::{Entity, ExperimentalMethod, is_cif_file, is_pdb_file, PDBError, Structure, UnitCell};

/// Holds all the data describing a macromolecular deposit, parsed from either an mmCIF or PDB file.
pub struct Deposit {
    /// classifies the molecule(s)
    ///
    /// This field should contain one of classifications from a curated list available at the [wwPDB website](http://www.wwpdb.org/)
    pub classification: Option<String>,
    /// deposition date
    pub dep_date: Option<String>,
    /// placeholder for keywords, which may be empty
    pub keywords: Vec<String>,
    /// four-character PDB code of this deposit, such as `2GB1` or `4HHB`
    pub id_code: String,
    /// title for a PDB entry
    ///
    /// This value is extracted either from a `TITLE` record of a PDB-formatted file
    /// or from a "_struct.title" entry of an mmCIF data.
    ///
    /// See  the [official documentation of the `TITLE` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#TITLE) for details
    pub title: Option<String>,
    /// describes how this structure was determined experimentally
    pub methods: Vec<ExperimentalMethod>,
    /// experimental resolution, when available
    pub resolution: Option<f64>,
    /// R-factor value, when available
    pub r_factor: Option<f64>,
    /// R-free value, when available
    pub r_free: Option<f64>,
    /// unit cell parameters, when available
    pub unit_cell: Option<UnitCell>,
    pub(crate) entities: HashMap<String, Entity>,
    pub(crate) structure: Structure,
}

impl Deposit {

    /// Creates a new, empty deposit for a given ``id_code``
    pub fn new(id_code: &str) -> Self {
        Deposit{
            classification: None,
            dep_date: None,
            keywords: vec![],
            id_code: id_code.to_string(),
            title: None,
            methods: vec![],
            resolution: None,
            r_factor: None,
            r_free: None,
            unit_cell: None,
            entities: Default::default(),
            structure: Structure::new(id_code),
        }
    }

    /// Detects the file format and parses its content into a [`Deposit`](Deposit)  struct.
    ///
    /// The method can recognise either mmCIF or PDB file format.
    pub fn from_file(file_name: &str) -> Result<Deposit, PDBError> {
        if is_cif_file(file_name)? { return Deposit::from_cif_file(file_name); }
        if is_pdb_file(file_name)? { return Deposit::from_pdb_file(file_name); }
        return Err(PDBError::InvalidFileFormat { file_name: file_name.to_string() });
    }

    /// Provides the number of entities of this [`Deposit`](Deposit)
    pub fn count_entities(&self) -> usize { self.entities.len() }

    /// Provides an iterator over the `entities` map.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_pdb::{EntityType, PDBError, Deposit};
    /// # fn main() -> Result<(), PDBError> {
    /// use bioshell_pdb::EntityType::Polymer;
    /// use bioshell_pdb::PolymerEntityType::{DNA, PolypeptideL};
    /// let cif_data = include_str!("../tests/test_files/5edw.cif");
    /// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
    /// // --- 5EDW deposit has five entities, including three polymer entities: two DNA chains and a single protein chain
    /// let mut polymer_entities = 0;
    /// for (id, entity) in deposit.entities() {
    ///     if entity.entity_type() == Polymer(DNA) || entity.entity_type() == Polymer(PolypeptideL) {
    ///         polymer_entities += 1;
    ///     }
    /// }
    /// # assert_eq!(polymer_entities, 3);
    /// # Ok(())
    /// # }
    /// ```
    ///
    pub fn entities(&self) -> impl Iterator<Item = (&String, &Entity)> { self.entities.iter() }

    /// Provides information about a given entity
    pub fn entity(&self, entity_id: &str) -> &Entity { &self.entities[entity_id] }

    /// returns a [`Structure`] object
    pub fn structure(&self) -> Structure { self.structure.clone() }

    /// Count models of a macromolecular structure(s) stored in this deposit
    pub fn count_models(&self) -> usize { self.structure.count_models() }
}