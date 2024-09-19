use std::collections::HashMap;
use std::io;
use std::io::{BufRead};
use std::time::Instant;
use log::{debug, info, warn};
use bioshell_cif::{cif_columns_by_name, CifLoop, read_cif_buffer, CifError, parse_item_or_error, value_or_default, entry_has_value, CifData};
use crate::{ExperimentalMethod, PdbAtom, PDBError, Structure, UnitCell, value_or_missing_key_pdb_error};
use crate::calc::Vec3;
use bioshell_cif::CifError::{ExtraDataBlock, MissingCifLoopKey, ItemParsingError, MissingCifDataKey};
use bioshell_io::open_file;
use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager, StandardResidueType};
use crate::crate_utils::find_deposit_file_name;
use crate::PDBError::{CifParsingError, IncorrectCompoundTypeName};

pub fn load_cif_reader<R: BufRead>(reader: R) -> Result<Structure, PDBError> {
    let start = Instant::now();
    let mut pdb_structure = Structure::new("");

    let cif_data = read_cif_buffer(reader)?;
    if cif_data.len() > 1 {
        return Err(CifParsingError(ExtraDataBlock));
    }
    let cif_data_block = &cif_data[0];

    if let Some(atoms_loop) = cif_data_block.first_loop("_atom_site.id") {
        cif_columns_by_name!(PdbAtomData, "_atom_site.id",
            "_atom_site.label_atom_id", "_atom_site.label_alt_id",
            "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id",
            "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.type_symbol",
            "_atom_site.pdbx_PDB_model_num", "_atom_site.auth_seq_id",
        );

        let extractor = PdbAtomData::new(atoms_loop)?;      // extracts atom data from each loop entry
        let mut tokens = [""; 15];                              // stores tokens from each extraction
        let mut model_ids: HashMap<usize, usize> = HashMap::new();       // links model id to model index in the model_coordinates structure
        for row in atoms_loop.rows() {
            extractor.data_items(&row, &mut tokens);
            let model_id = tokens[13].parse::<usize>().unwrap();  // model_id of the current atom
            if ! model_ids.contains_key(&model_id) {                    // if we have a new model...
                model_ids.insert(model_id, model_ids.len());
                pdb_structure.model_coordinates.push(vec![]);
            }
            let mdl_idx = model_ids[&model_id];
            let x = tokens[7].parse::<f64>().unwrap();
            let y = tokens[8].parse::<f64>().unwrap();
            let z = tokens[9].parse::<f64>().unwrap();
            let pos = Vec3::new(x, y, z);

            if model_ids.len() == 1 {
                pdb_structure.model_coordinates[mdl_idx].push(pos.clone());
                match create_pdb_atom(&tokens, pos) {
                    Ok(a) => { pdb_structure.atoms.push(a); }
                    Err(e) => {
                        let data_line = row.join(" ");
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
                pdb_structure.model_coordinates[mdl_idx].push(pos);
            }
        }
    } else {
        warn!("mmCIF data has no atoms: double check if it contains _atom_site loop block");
        return Err(CifParsingError(MissingCifLoopKey{ item_key: "_atom_site".to_string() }));
    }

    // --- header data
    pdb_structure.classification = cif_data_block.get_item("_struct_keywords.pdbx_keywords");
    pdb_structure.id_code = value_or_missing_key_pdb_error!(cif_data_block, "_entry.id", String);
    pdb_structure.dep_date = cif_data_block.get_item("_pdbx_database_status.recvd_initial_deposition_date");
    pdb_structure.title = cif_data_block.get_item("_struct.title");

    // --- exp details and resolution
    pdb_structure.methods = ExperimentalMethod::from_cif_data(cif_data_block);
    pdb_structure.resolution = cif_data_block.get_item("_refine.ls_d_res_high");
    pdb_structure.r_factor = cif_data_block.get_item("_refine.ls_R_factor_obs");
    pdb_structure.r_free = cif_data_block.get_item("_refine.ls_R_factor_R_free");

    // --- crystallography parameters
    pdb_structure.unit_cell = if let Ok(uc) = UnitCell::from_cif_data(cif_data_block) { Some(uc) } else { None };
    pdb_structure.update();
    debug!("Structure loaded in: {:?}", start.elapsed());

    return Ok(pdb_structure);
}

/// Reads a [`Structure`](Structure) from a mmCIF file
///
pub fn load_cif_file(file_name: &str) -> Result<Structure, PDBError> {

    info!("Loading an mmCIF deposit: {}", file_name);
    let reader = open_file(file_name)?;

    return load_cif_reader(reader);
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

/// A helper function to create an atom based on given string tokens
fn create_pdb_atom(tokens: &[&str; 15], pos: Vec3) -> Result<PdbAtom, CifError> {

    let serial = parse_item_or_error!(tokens[0], i32);
    let name = tokens[1].to_string();
    let alt_loc = value_or_default(tokens[2], ' ');
    let res_name = tokens[3].to_string();
    let chain_id = tokens[4].to_string();
    let res_seq = if entry_has_value(tokens[5]) {
        parse_item_or_error!(tokens[5], i32)
    } else {
        parse_item_or_error!(tokens[14], i32)
    };

    let i_code = value_or_default(tokens[6], ' ');
    let occupancy = parse_item_or_error!(tokens[10], f64);
    let temp_factor = parse_item_or_error!(tokens[11], f64);
    let element = Some(tokens[12].to_string());
    let charge = None;
    let is_hetero_atom = false;
    let secondary_struct_type: u8 = 12;
    let a = PdbAtom {
        serial, name, alt_loc, res_name,
        chain_id, res_seq, i_code, pos, occupancy, temp_factor,
        element, charge, is_hetero_atom, secondary_struct_type,
    };

    return Ok(a);
}

fn load_residue_types(cif_data_block: &CifData) -> Result<(), PDBError>{
    if let Some(monomers_loop) = cif_data_block.first_loop("_chem_comp.id") {
        cif_columns_by_name!(MonomersData, "_chem_comp.id", "_chem_comp.type",);

        let extractor = MonomersData::new(monomers_loop)?;
        let mut tokens = [""; 2];
        let mut rts = ResidueTypeManager::get();
        for row in monomers_loop.rows() {
            extractor.data_items(&row, &mut tokens);
            match  MonomerType::try_from(tokens[1]) {
                Ok(chem_type) => {
                    let rt = ResidueType::from_attrs(tokens[0], StandardResidueType::UNK, chem_type);
                    rts.register_residue_type(rt);
                }
                Err(_) => {return Err(IncorrectCompoundTypeName{ compound_id: tokens[0].to_string(), compound_type: tokens[1].to_string()})}
            };
        }
    }
    Ok(())
}

