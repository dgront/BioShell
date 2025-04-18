use std::collections::{HashMap, HashSet};
use std::path::Path;
use bioshell_cif::CifData;
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
    /// the number of structural models in this deposit
    pub n_models: usize,
    pub(crate) entities: HashMap<String, Entity>,
    pub(crate) structure: Option<Structure>,
    pub(crate) cif_buffer: Option<CifData>
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
            n_models: 0,
            entities: Default::default(),
            structure: None,
            cif_buffer: None,
        }
    }

    /// Detects the file format and parses its content into a [`Deposit`](Deposit)  struct.
    ///
    /// The method can recognise either mmCIF or PDB file format.
    ///
    /// # Examples
    /// ```
    /// # use bioshell_pdb::Deposit;
    /// # use bioshell_pdb::PDBError;
    /// # fn main() -> Result<(), PDBError> {
    /// use std::path::{Path, PathBuf};
    /// let deposit = Deposit::from_file("tests/test_files/2gb1.cif")?;       // from &str
    /// let deposit = Deposit::from_file(String::from("tests/test_files/2gb1.cif"));  // from String
    /// let deposit = Deposit::from_file(Path::new("tests/test_files/2gb1.cif"));     // from &Path
    /// let deposit = Deposit::from_file(PathBuf::from("tests/test_files/2gb1.cif"));
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_file<P: AsRef<Path>>(file_path: P) -> Result<Deposit, PDBError> {

        if is_cif_file(&file_path)? { return Deposit::from_cif_file(&file_path); }
        if is_pdb_file(&file_path)? { return Deposit::from_pdb_file(&file_path); }
        return Err(PDBError::InvalidFileFormat { file_name: file_path.as_ref().to_string_lossy().into_owned() });
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

    /// Returns an iterator over all unique chain IDs present in the deposit.
    ///
    /// This method collects all chain identifiers (`chain_ids`) from each `Entity`
    /// stored in the deposit and returns an iterator over the unique identifiers.
    ///
    /// # Returns
    /// An iterator over `&str` representing the unique chain IDs.
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::Deposit;
    /// # use bioshell_pdb::PDBError;
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/5edw.cif");
    /// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
    /// let chain_ids: Vec<_> = deposit.chain_ids().collect();
    /// for id in &chain_ids {
    ///     println!("{}", id);
    /// }
    /// assert_eq!(chain_ids.len(), 3);
    /// # Ok(())
    /// # }
    /// ```
    pub fn chain_ids(&self) -> impl Iterator<Item = &str> {

        let mut unique_ids = HashSet::new();
        for entity in self.entities.values() {
            for id in entity.chain_ids() {
                unique_ids.insert(id.as_str());
            }
        }
        unique_ids.into_iter()
    }

    /// Provides information about a given entity
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::Deposit;
    /// # use bioshell_pdb::PDBError;
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/2gb1.cif");
    /// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
    /// let entity = deposit.entity("1");
    /// assert!(entity.is_some());
    /// let entity = deposit.entity("2");   // 2gb1 deposit has only one entity!
    /// assert!(entity.is_none());
    /// # Ok(())
    /// # }
    /// ```
    pub fn entity(&self, entity_id: &str) -> Option<&Entity>  { self.entities.get(entity_id) }

    /// Returns a [`Structure`] object.
    ///
    /// A structure is lazily parsed from a PDB or mmCIF file, i.e. it is not parsed until this method is called.
    /// This allows fast access to a deposit's data, without the need to parse the whole file.
    /// If the structure has been already parsed, the method returns clone of a [`Structure`] object.
    pub fn structure(&self) -> Option<Structure> {
        match &self.structure {
            Some(structure) => Some(structure.clone()),
            None => {
                match &self.cif_buffer {
                    Some(cif_data) => Self::structure_from_cif_data(cif_data).ok(),
                    None => None,
                }
            },
        }
    }
}