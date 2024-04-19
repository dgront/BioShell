use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use log::{debug, info};
use crate::{ExperimentalMethod, PdbAtom, PdbHeader, PdbHelix, PDBRemarks, PdbSheet, PdbTitle, residue_id_from_ter_record, ResidueId, SecondaryStructureTypes, Structure};
use crate::pdb_atom_filters::{ByResidueRange, PdbAtomPredicate};
use crate::pdb_parsing_error::PDBError;

/// Reads PDB-formatted content from a buffer.
///
/// This function allows reading PDB structures from `String`
/// # Example
/// ```
/// use bioshell_pdb::load_pdb_reader;
/// use std::io::BufReader;
/// let pdb_txt: &str =
/// "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
/// ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C
/// ATOM     35  CA  TYR A   3      -7.173  -2.314   2.811  1.00  0.08           C
/// ATOM     56  CA  LYS A   4      -3.922  -3.881   4.044  1.00  0.10           C
/// ATOM     78  CA  LEU A   5      -0.651  -2.752   2.466  1.00  0.11           C
/// ATOM     97  CA  ILE A   6       2.338  -5.105   2.255  1.00  0.13           C";
///
/// let strctr = load_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
/// let seq = strctr.sequence("A");
/// assert_eq!(seq.to_string(), "MTYKLI");
/// ```
pub fn load_pdb_reader<R: BufRead>(reader: R) -> Result<Structure, PDBError> {

    let start = Instant::now();
    let mut pdb_structure = Structure::new();

    let mut atoms: Vec<PdbAtom> = vec![];
    let mut helices: Vec<PdbHelix> = vec![];
    let mut strands: Vec<PdbSheet> = vec![];
    let mut remarks = PDBRemarks::new();

    for line in reader.lines() {
        let line = line?;

        if line.len() < 6 { continue }  // --- remove empty lines or those without a valid record
        // Check that the line has a valid PDB record type
        let record = line[0..6].trim();
        match record {
            "TER" => {
                if let Ok(ter_res) = residue_id_from_ter_record(&line) {
                    let ter_chain = ter_res.chain_id.clone();
                    pdb_structure.ter_atoms.insert(ter_chain, ter_res);
                } else {                    // --- assign that TER to the very last atom of the current chain
                    if let Some(last_atom) = atoms.last() {
                        pdb_structure.ter_atoms.insert(last_atom.chain_id.clone(), ResidueId::try_from(last_atom).unwrap());
                    }
                }
            }
            "HEADER" => {
                let header = PdbHeader::new(&line);
                pdb_structure.header = Some(header);
            },
            "EXPDTA" => {
                let header = PdbHeader::new(&line);
                pdb_structure.methods = ExperimentalMethod::from_expdata_line(&line);
            },
            "TITLE" => {
                if pdb_structure.title == None {
                    pdb_structure.title = Some(PdbTitle::new(&line));
                } else {
                    let title = pdb_structure.title.as_mut().unwrap();
                    title.append(&line);
                }
            },
            "HELIX" => {
                let helix = PdbHelix::from_helix_line(&line);
                helices.push(helix);
            }
            "SHEET" => {
                let strand = PdbSheet::from_sheet_line(&line);
                strands.push(strand);
            }
            "ATOM" => {
                atoms.push(PdbAtom::from_atom_line(&line));
            },
            "HETATM" => {
                atoms.push(PdbAtom::from_atom_line(&line));
            },
            "ENDMDL" => {
                break;
            },
            "REMARK" => {
                remarks.add_remark(&line);
            }
            _ => {},
        };
    }
    debug!("{:} atoms loaded",atoms.len());

    // ---------- Annotate secondary structure
    for i in 0..helices.len() {
        let from = helices[i].init_res_id();
        let to = helices[i].end_res_id();
        let check = ByResidueRange::new(from,to);
        for a in &mut atoms {
            if check.check(a) { a.secondary_struct_type = helices[i].helix_class }
        }
    }
    for i in 0..strands.len() {
        let from = strands[i].init_res_id();
        let to = strands[i].end_res_id();
        let check = ByResidueRange::new(from,to);
        for a in &mut atoms {
            if check.check(a) { a.secondary_struct_type = SecondaryStructureTypes::Strand as u8 }
        }
    }

    // ---------- Extract values stored in remarks
    pdb_structure.resolution = remarks.resolution();

    pdb_structure.atoms = atoms;
    pdb_structure.update();
    debug!("Structure loaded in: {:?}", start.elapsed());

    Ok(pdb_structure)
}

/// Reads a [`Structure`](Structure) from a PDB file
///
pub fn load_pdb_file(file_name: &str) -> Result<Structure, PDBError> {

    info!("Loading a PDB deposit: {}", file_name);

    let file = File::open(file_name)?;
    let reader = BufReader::new(file);
    return load_pdb_reader(reader);
}