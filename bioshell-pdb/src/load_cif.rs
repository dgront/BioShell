use std::io;
use std::io::{BufRead};
use std::time::Instant;
use log::{debug, info};
use bioshell_cif::{read_cif_buffer, CifData, CifTable};
use crate::{Deposit, Entity, ExperimentalMethod, PdbAtom, PDBError, SecondaryStructureTypes, Structure, UnitCell};
use bioshell_cif::CifError::{ExtraDataBlock};
use bioshell_io::open_file;
use crate::pdb_sheet::PdbSheet;
use crate::pdb_helix::PdbHelix;
use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager, StandardResidueType};
use crate::crate_utils::find_deposit_file_name;
use crate::ExperimentalMethod::{ElectronCrystallography, FiberDiffraction, ElectronMicroscopy, XRay};
use crate::pdb_atom_filters::{ByResidueRange, PdbAtomPredicate};
use crate::PDBError::{CifParsingError, IncorrectCompoundTypeName};

impl Deposit {

    /// Reads a [`Deposit`](Deposit) from a mmCIF data buffer.
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_pdb::{Deposit, PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/2gb1.cif");
    /// let reader = BufReader::new(cif_data.as_bytes());
    /// let deposit = Deposit::from_cif_reader(reader)?;
    /// assert_eq!(deposit.count_entities(), 1);
    /// Ok(())
    /// # }
    /// ```
    pub fn from_cif_reader<R: BufRead>(reader: R) -> Result<Deposit, PDBError> {
        let start = Instant::now();
        // --- parse the file content into CIF struct
        let mut cif_data = read_cif_buffer(reader)?;
        if cif_data.len() > 1 {
            return Err(CifParsingError(ExtraDataBlock));
        }
        let cif_data_block = &cif_data[0];

        // --- name of the data block must be the deposit ID
        let mut deposit = Deposit::new(cif_data_block.name());

        // ---------- load the residue types into the ResidueTypeManager before atoms and entities
        load_residue_types(cif_data_block)?;

        // --- header data
        deposit.classification = cif_data_block.get_item("_struct_keywords.pdbx_keywords");
        deposit.dep_date = cif_data_block.get_item("_pdbx_database_status.recvd_initial_deposition_date");
        deposit.title = cif_data_block.get_item("_struct.title");

        // --- exp details and resolution
        deposit.methods = ExperimentalMethod::from_cif_data(cif_data_block);
        if deposit.methods.contains(&XRay) {
            deposit.resolution = cif_data_block.get_item("_refine.ls_d_res_high");
            deposit.r_factor = cif_data_block.get_item("_refine.ls_R_factor_obs");
            deposit.r_free = cif_data_block.get_item("_refine.ls_R_factor_R_free");
        }
        if deposit.methods.contains(&ElectronMicroscopy) {
            deposit.resolution = cif_data_block.get_item("_em_3d_reconstruction.resolution");
        }
        if deposit.methods.contains(&FiberDiffraction) {
            deposit.resolution = cif_data_block.get_item("_pd_proc.reflns_resolution");
        }
        if deposit.methods.contains(&ElectronCrystallography) {
            deposit.resolution = cif_data_block.get_item("_em_3d_reconstruction.resolution");
        }
        if let Some(keywds) = cif_data_block.get_item::<String>("_struct_keywords.text") {
            deposit.keywords = keywds.split(",").map(|s| s.to_string()).collect();
        }

        // --- entity stuff
        for (entity_id, entity) in Entity::from_cif_data(cif_data_block).unwrap() {
            deposit.entities.insert(entity_id, entity);
        }

        // --- crystallography parameters
        deposit.unit_cell = if let Ok(uc) = UnitCell::from_cif_data(cif_data_block) { Some(uc) } else { None };
        debug!("{} deposit loaded in: {:?}", &deposit.id_code, start.elapsed());

        // --- store the CIF data block for future access
        deposit.cif_buffer = Some(cif_data.remove(0));

        return Ok(deposit);
    }

    /// Reads a [`Deposit`](Deposit) from a mmCIF file.
    ///
    /// This method loads all the data from the CIF file, except for the atoms; these are loaded
    /// in a lazy manner at the first at the first time the [`Structure`](Structure) is accessed.
    ///
    pub fn from_cif_file(file_name: &str) -> Result<Deposit, PDBError> {
        info!("Loading an mmCIF deposit: {}", file_name);
        let reader = open_file(file_name)?;

        return Self::from_cif_reader(reader);
    }

    /// Creates a new [`Structure`](Structure) from a CifData block stored in this [`Deposit`](Deposit).
    ///
    /// This method is called by the [`Deposit::structure()`](Deposit::structure()) method.
    pub(crate) fn structure_from_cif_data(cif_data_block: &CifData) -> Result<Structure, PDBError> {

        let start = Instant::now();

        // --- Create a new structure
        let mut structure = Structure::new(cif_data_block.name());

        // --- load all the atoms
        PdbAtom::from_cif_data(cif_data_block, &mut structure)?;


        // todo: fix the helix type! Now it's only an alpha helix
        // --- annotate secondary structure
        let helices = PdbHelix::from_cif_data(cif_data_block)?;
        for h in &helices {
            let range = ByResidueRange::new(h.init_res_id(), h.end_res_id());
            structure.atoms.iter_mut().for_each(|a| if range.check(a) {
                a.secondary_struct_type = SecondaryStructureTypes::RightAlphaHelix
            });
        }
        let strands = PdbSheet::from_cif_data(cif_data_block)?;
        for s in &strands {
            let range = ByResidueRange::new(s.init_res_id(), s.end_res_id());
            structure.atoms.iter_mut().for_each(|a| if range.check(a) {
                a.secondary_struct_type = SecondaryStructureTypes::Strand
            });
        }
        structure.update();
        debug!("Structure loaded in: {:?}", start.elapsed());

        return Ok(structure);
    }
}

/// Returns true if a given file is in CIF format.
///
/// This function simply tests whether the first non-empty data line of a given file starts with ``data_``,
/// Otherwise, it returns ``false``. When the file can't be open returns I/O error.
///
/// # Examples
/// ```
/// use bioshell_cif::is_cif_file;
/// let try_2gb1 = is_cif_file("./tests/test_files/2gb1.cif");
/// assert!(try_2gb1.is_ok());
/// assert!(try_2gb1.unwrap());
/// let try_2gb1 = is_cif_file("./tests/test_files/2gb1.pdb");
/// assert!(try_2gb1.is_ok());
/// assert!(!try_2gb1.unwrap());
/// ```
pub fn is_cif_file(file_path: &str) -> io::Result<bool> {
    let reader = open_file(file_path)?;

    let cif_starts_with = ["data_"];
    for line in reader.lines() {
        let line = line?;
        if !line.is_empty() {
            return Ok(cif_starts_with.iter().any(|s|line.starts_with(s)));
        }
    }

    return Ok(false);
}

static CIF_PREFIXES: [&str; 2] = ["", "pdb"];
static CIF_SUFFIXES: [&str; 6] = [".cif", ".cif.gz", ".gz", ".CIF", ".CIF.gz", ""];

/// Attempts to find a CIF file in a given directory.
///
/// Looks in the specified path for a file with a given PDB data, identified by
/// a given PDB code. For a given 4-character ID (digit + 3 letters), the method checks
/// the following possibilities:
///
/// - `given_path/1abc`
/// - `given_path/1ABC`
/// - `given_path/1abc.cif`
/// - `given_path/1ABC.cif`
/// - `given_path/1ABC.CIF`
/// - `given_path/pdb1abc`
/// - `given_path/PDB1ABC`
/// - `given_path/pdb1abc.cif`
/// - `given_path/pdb1abc.cif.gz`
/// - `given_path/ab/pdb1abc.cif`
/// - `given_path/ab/pdb1abc.cif.gz`
///
/// where `1abc` and `1ABC` denote a lower-case and an upper-case PDB ID, respectively. Returns
/// the name of the PDB file that was found or an error.
///
/// # Arguments
///
/// * `pdb_code` - A four-character PDB ID.
/// * `pdb_path` - Directory to look into.
///
/// # Example
/// ```
/// use bioshell_pdb::find_cif_file_name;
/// let result = find_cif_file_name("2gb1", "./tests/test_files/");
/// assert!(result.is_ok());
/// assert_eq!(result.unwrap(), "./tests/test_files/2gb1.cif");
/// ```
///
pub fn find_cif_file_name(pdb_code: &str, pdb_path: &str) -> Result<String, io::Error> {
    find_deposit_file_name(pdb_code, pdb_path, &CIF_PREFIXES, &CIF_SUFFIXES)
}



/// Loads residue types from a CIF file directly into the [`ResidueTypeManager`]
fn load_residue_types(cif_data_block: &CifData) -> Result<(), PDBError>{
    let monomers = CifTable::new(cif_data_block, "_chem_comp.",["id", "type",])?;
    let mut rts = ResidueTypeManager::get();
    for tokens in monomers.iter() {
        match  MonomerType::try_from(tokens[1]) {
            Ok(chem_type) => {
                let rt = ResidueType::from_attrs(tokens[0], StandardResidueType::UNK, chem_type);
                rts.register_residue_type(rt);
            }
            Err(_) => {return Err(IncorrectCompoundTypeName{ compound_id: tokens[0].to_string(), compound_type: tokens[1].to_string()})}
        };
    }
    Ok(())
}


