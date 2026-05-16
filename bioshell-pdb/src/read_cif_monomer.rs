use std::io::BufRead;
use bioshell_cif::{CifError, CifTable, read_cif_buffer};
use crate::{format_atom_name, PdbAtom, PDBError, Structure};
use crate::calc::Vec3;

/// Reads a library of small molecules (e.g. ligands) into a  ([`Structure`](Structure)
///
/// The function reads each CIF data block assuming it contains a separate small molecule and places it
/// into a separate residue in the chain `'A'` of the structure.
///
/// # Example
/// ```
/// let ala = include_str!("../tests/test_files/ala.cif");
/// let reader = std::io::BufReader::new(ala.as_bytes());
/// let strctr = bioshell_pdb::read_cif_monomers(reader, None).unwrap();
/// assert_eq!(strctr.count_atoms(), 13);
/// ```
pub fn read_cif_monomers<R: BufRead>(reader: R, molecule_name: Option<&str>) -> Result<Structure, PDBError> {

    let blocks = read_cif_buffer(reader)?;
    if blocks.len() < 1 { return Err(PDBError::NoStructureDataLoaded); }
    let molecule_name = molecule_name.unwrap_or(blocks[0].name());
    let mut strctr = Structure::new(molecule_name);
    for (b_id, b) in blocks.iter().enumerate() {
        let atom_table = CifTable::new(b, "_chem_comp_atom",
            ["pdbx_ordinal","atom_id", "pdbx_component_comp_id", "model_Cartn_x", "model_Cartn_y", "model_Cartn_z", "charge", "type_symbol"])?;
        for [idx, name, res_name, x, y, z, charge, elem] in atom_table.iter() {
            let mut a = PdbAtom::new();
            a.serial = idx.parse::<i32>().map_err(|_e| CifError::CantParseIntValue { value: idx.to_string() })?;
            a.entity_id = "A".to_string();
            a.chain_id = "A".to_string();
            a.res_seq = b_id as i32 + 1;
            a.res_name = res_name.to_string();
            a.name = format_atom_name(name, Some(elem));
            if charge!= "0" {
                a.charge = Some(charge.to_string());
            }
            a.element = Some(elem.to_string());
            let vx = x.parse::<f64>().map_err(|_e| CifError::CantParseFloatValue { value: x.to_string() })?;
            let vy = y.parse::<f64>().map_err(|_e| CifError::CantParseFloatValue { value: x.to_string() })?;
            let vz = z.parse::<f64>().map_err(|_e| CifError::CantParseFloatValue { value: x.to_string() })?;
            a.pos = Vec3::new(vx, vy, vz);
            strctr.push_atom(a);
        }
    }

    return Ok(strctr);
}
