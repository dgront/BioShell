use std::fs::File;
use std::io::{BufRead, BufReader};
use std::str::FromStr;
use std::time::Instant;
use log::{debug, info, warn};
use bioshell_cif::{cif_columns_by_name, CifLoop, read_cif_buffer, CifError};
use crate::{ExperimentalMethod, PdbAtom, PdbHeader, Structure};
use crate::calc::Vec3;
use bioshell_cif::CifError::{ExtraDataBlock, MissingCifDataKey, MissingCifLoopKey};

pub fn load_cif_reader<R: BufRead>(reader: R) -> Result<Structure, CifError> {
    let start = Instant::now();
    let mut pdb_structure = Structure::new();

    let cif_data = read_cif_buffer(reader);
    if cif_data.len() > 1 {
        return Err(ExtraDataBlock);
    }
    let cif_data_block = &cif_data[0];

    if let Some(atoms_loop) = cif_data_block.first_loop("_atom_site") {
        cif_columns_by_name!(PdbAtomData, "_atom_site.id",
            "_atom_site.label_atom_id", "_atom_site.label_alt_id",
            "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_seq_id",
            "_atom_site.pdbx_PDB_ins_code",
            "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
            "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.type_symbol",
            "_atom_site.pdbx_PDB_model_num",
        );

        let extractor = PdbAtomData::new(atoms_loop)?;
        let mut tokens = [""; 14];
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

                let serial = tokens[0].parse::<i32>().unwrap();
                let name = tokens[1].to_string();
                let alt_loc = value_or_default(tokens[2], ' ');
                let res_name = tokens[3].to_string();
                let chain_id = tokens[4].to_string();
                let res_seq = tokens[5].parse::<i32>().unwrap();
                let i_code = value_or_default(tokens[6], ' ');
                let occupancy = tokens[10].parse::<f64>().unwrap();
                let temp_factor = tokens[11].parse::<f64>().unwrap();
                let element = Some(tokens[12].to_string());
                let charge = None;
                let is_hetero_atom = false;
                let secondary_struct_type: u8 = 12;
                let a = PdbAtom{
                    serial, name, alt_loc, res_name,
                    chain_id, res_seq, i_code, pos, occupancy, temp_factor,
                    element, charge, is_hetero_atom, secondary_struct_type,
                };
                pdb_structure.atoms.push(a);
            } else {
                pdb_structure.model_coordinates[model_id-1].push(pos);
            }
        }
    } else {
        warn!("mmCIF data has no atoms: double check if it contains _atom_site loop block");
        return Err(MissingCifLoopKey{ item_key: "_atom_site".to_string() });
    }


    // --- header data
    pdb_structure.header = PdbHeader::from_cif(cif_data_block);
    let methods: Option<String> = cif_data_block.get_item("_exptl.method");
    match methods {
        Some(methods) => pdb_structure.methods = ExperimentalMethod::from_expdata_line(&methods),
        None => return Err(MissingCifDataKey{ item_key: "_exptl.method".to_string() })
    }
    // --- resolution
    pdb_structure.resolution = cif_data_block.get_item("_refine.ls_d_res_high");
    pdb_structure.update();
    debug!("Structure loaded in: {:?}", start.elapsed());

    return Ok(pdb_structure);
}

/// Reads a [`Structure`](Structure) from a PDB file
///
pub fn load_cif_file(file_name: &str) -> Result<Structure, CifError> {

    info!("Loading an mmCIF deposit: {}", file_name);

    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    return load_cif_reader(reader);
}



fn value_or_default<T: FromStr>(data_entry: &str, default_val: T) -> T {
    if data_entry == "?" || data_entry == "." { return default_val }
    return data_entry.parse().ok().unwrap();
}

