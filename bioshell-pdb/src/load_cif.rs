use std::io::{BufRead};
use std::time::Instant;
use log::{debug, info, warn};
use bioshell_cif::{cif_columns_by_name, CifLoop, read_cif_buffer, CifError, parse_item_or_error, value_or_default, entry_has_value};
use crate::{ExperimentalMethod, PdbAtom, PDBError, Structure, UnitCell, value_or_missing_key_pdb_error};
use crate::calc::Vec3;
use bioshell_cif::CifError::{ExtraDataBlock, MissingCifLoopKey, ItemParsingError, MissingCifDataKey};
use bioshell_io::open_file;
use crate::PDBError::CifParsingError;

pub fn load_cif_reader<R: BufRead>(reader: R) -> Result<Structure, PDBError> {
    let start = Instant::now();
    let mut pdb_structure = Structure::new("");

    let cif_data = read_cif_buffer(reader);
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

        let extractor = PdbAtomData::new(atoms_loop)?;
        let mut tokens = [""; 15];
        for row in atoms_loop.rows() {
            extractor.data_items(&row, &mut tokens);
            let model_id = tokens[13].parse::<usize>().unwrap();
            let x = tokens[7].parse::<f64>().unwrap();
            let y = tokens[8].parse::<f64>().unwrap();
            let z = tokens[9].parse::<f64>().unwrap();
            let pos = Vec3::new(x, y, z);

            if pdb_structure.model_coordinates.len() == model_id - 1 {
                pdb_structure.model_coordinates.push(vec![]);
            }

            if model_id == 1 {
                pdb_structure.model_coordinates[model_id-1].push(pos.clone());
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
                pdb_structure.model_coordinates[model_id-1].push(pos);
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
