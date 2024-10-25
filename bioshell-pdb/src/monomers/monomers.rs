use std::collections::HashMap;
use std::num::ParseIntError;
use std::sync::{Mutex, MutexGuard};
use clap::__macro_refs::once_cell::sync::Lazy;
use log::{debug, info};
use bioshell_cif::{CifError, CifTable, parse_bool, read_cif_file};
use bioshell_io::find_bioshell_path;
use bioshell_seq::chemical::{ResidueType, ResidueTypeManager, ResidueTypeProperties, StandardResidueType};
use bioshell_seq::chemical::StandardResidueType::{UNK, UNL, GAP, GPE};

/// Defines an atom of a monomer structure
pub struct MonomerAtom {
    /// The index of the atom in the monomer
    pub atom_number: i32,
    /// The name of the atom
    pub atom_name: String,
    /// Leaving atoms are removed from a monomer while a polymer is formed
    pub is_leaving: bool,
    /// Chemical element of the atom
    pub atom_type: String,
}

/// Defines a monomer chemical structure
pub struct Monomer {
    residue_type: ResidueType,
    atoms: Vec<MonomerAtom>,
}

impl Monomer {
    /// Provides the reference to the definition of this residue type
    pub fn residue_type(&self) -> &ResidueType { &self.residue_type }

    /// Provides the reference to the list of atoms in this residue
    pub fn atoms(&self) -> &Vec<MonomerAtom> { &self.atoms }

    /// Total number of atoms in this monomer, including leaving atoms
    ///
    /// A leaving atom is an atom that is present in the structure of the monomer, but is removed
    /// when a polymer is built, e.g. "OXT" in the case of amino acids.
    ///
    /// ```
    /// use bioshell_pdb::monomers::MonomerManager;
    /// let manager = MonomerManager::get();
    /// let ala = manager.by_code3("ALA").unwrap();
    /// assert_eq!(ala.count_all_atoms(), 13);
    /// ```
    pub fn count_all_atoms(&self) -> usize { self.atoms.len() }

    /// The number of atoms in this monomer residue, excluding leaving atoms
    ///
    /// ```
    /// use bioshell_pdb::monomers::MonomerManager;
    /// let manager = MonomerManager::get();
    /// let ala = manager.by_code3("ALA").unwrap();
    /// assert_eq!(ala.count_residue_atoms(), 10);
    /// ```
    pub fn count_residue_atoms(&self) -> usize { self.atoms.iter().filter(|a| !a.is_leaving).count() }


    /// The number of heavy atoms in this monomer residue, excluding leaving atoms and hydrogens
    ///
    /// ```
    /// use bioshell_pdb::monomers::MonomerManager;
    /// let manager = MonomerManager::get();
    /// let ala = manager.by_code3("ALA").unwrap();
    /// assert_eq!(ala.count_residue_heavy(), 5);
    /// ```
    pub fn count_residue_heavy(&self) -> usize {
        self.atoms.iter().filter(|a| ! (a.is_leaving ||  &a.atom_type == "H")).count()
    }

}

/// Provides a definition of a monomer type.
///
/// A respective structure must be registered in this manager by a user; respective files in CIF format
/// may be downloaded from the [PDB Ligand Repository](https://files.rcsb.org/ligands/), e.g.
/// here is the direct link to [alanine definition](https://files.rcsb.org/ligands/view/ALA.cif)
///
/// Standard monomers: 20 amino acids and 9 nucleotide variants are provided by default.
pub struct MonomerManager {
    by_code_3: HashMap<String, Monomer>
}

impl MonomerManager {

    pub(crate) fn new() -> MonomerManager {
        let mut mgr = MonomerManager { by_code_3: HashMap::new() };

        if let Some(path) = find_bioshell_path() {
            let path = path.join("bioshell-pdb").join("data").join("monomers");
            for rt in &StandardResidueType::TYPES {
                match rt {
                    UNK | UNL | GAP | GPE=> {}
                    _ => {
                        let fname = format!("{}.cif", rt.code3());
                        let out = mgr.load_cif_file(path.join(fname).to_str().unwrap());
                        if out.is_err() {
                            panic!("Failed to load residue structure for '{}'", rt.code3());
                        } else {
                            debug!("Loaded monomer '{}'", rt.code3());
                        }
                    }
                }
            }
        } else { panic!("Can't pre-load standard monomers!"); }

        info!("Loaded {} standard monomer residue structures", mgr.count());

        return mgr;
    }

    /// Provides a singleton instance of a [`MonomerManager`] object
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::monomers::MonomerManager;
    /// let rs_mgr = MonomerManager::get();
    /// assert!(rs_mgr.by_code3("ALA").is_some());
    /// assert_eq!(rs_mgr.count(), 29);
    /// ```
    pub fn get() -> MutexGuard<'static, MonomerManager> {
        KNOWN_MONOMERS.lock().unwrap()
    }

    /// Provides the residue structure for a given three-letter code.
    /// Returns an option that contains the monomer or None if it hasn’t been registered
    pub fn by_code3(&self, code_3: &str) -> Option<&Monomer> { self.by_code_3.get(code_3) }

    /// Counts the residue structures registered in this manager
    pub fn count(&self) -> usize { self.by_code_3.len() }

    pub fn load_cif_file(&mut self, file_name: &str) -> Result<(), CifError> {
        let cif_data = read_cif_file(file_name)?;
        for data_block in &cif_data {
            let res_name = data_block.name();
            let mut atoms: Vec<MonomerAtom> = vec![];
            let atom_table = CifTable::new(data_block, "_chem_comp_atom",
        ["pdbx_ordinal", "atom_id", "pdbx_leaving_atom_flag", "type_symbol"])?;
            for [idx, atom_name, is_leaving, element] in atom_table.iter() {
                let leaving_bool = parse_bool(is_leaving)?;
                atoms.push(MonomerAtom {
                    atom_number: idx.parse().map_err(|e: ParseIntError| CifError::ItemParsingError {
                        item: "pdbx_ordinal".to_string(),
                        type_name: "i32".to_string(),
                        details: e.to_string(),
                    })?,
                    atom_name: atom_name.to_string(),
                    is_leaving: leaving_bool,
                    atom_type: element.to_string(),
                });
            }
            if let Some(residue_type) = ResidueTypeManager::get().by_code3(res_name) {
                self.by_code_3.insert(res_name.to_string(), Monomer { residue_type: residue_type.clone(), atoms, });
            };
        }

        Ok(())
    }
}


pub static KNOWN_MONOMERS: Lazy<Mutex<MonomerManager>>
                = Lazy::new(|| Mutex::new(MonomerManager::new()));

