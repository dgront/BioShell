use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use log::{debug, info};
use bioshell_cif::read_cif_buffer;
use crate::{PDBError, PdbHeader, Structure};
use crate::PDBError::ExtraDataBlock;

pub fn load_cif_reader<R: BufRead>(reader: R) -> Result<Structure, PDBError> {
    let start = Instant::now();
    let mut pdb_structure = Structure::new();

    let mut cif_data = read_cif_buffer(reader);
    if cif_data.len() > 1 {
        return Err(ExtraDataBlock);
    }
    let cif_data_block = &cif_data[0];

    // --- header data
    pdb_structure.header = PdbHeader::from_cif(cif_data_block);
    // --- resolution
    pdb_structure.resolution = cif_data_block.get_item("_refine.ls_d_res_high");
    // pdb_structure.update();
    debug!("Structure loaded in: {:?}", start.elapsed());

    return Ok(pdb_structure);
}

/// Reads a [`Structure`](Structure) from a PDB file
///
pub fn load_cif_file(file_name: &str) -> Result<Structure, PDBError> {

    info!("Loading an mmCIF deposit: {}", file_name);

    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    return load_cif_reader(reader);
}