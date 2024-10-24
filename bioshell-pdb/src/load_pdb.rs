use std::collections::HashMap;
use std::io;
use std::io::{BufRead};
use std::time::Instant;
use log::{debug, info};
use bioshell_io::open_file;
use bioshell_seq::chemical::{ResidueTypeManager, ResidueTypeProperties};
use bioshell_seq::sequence::Sequence;
use crate::pdb_title::PdbTitle;
use crate::pdb_header::PdbHeader;
use crate::remarks::PDBRemarks;
use crate::{Deposit, ExperimentalMethod, PdbAtom, PdbHelix, PdbSheet, residue_id_from_ter_record, ResidueId, SecondaryStructureTypes, Structure, UnitCell};
use crate::calc::Vec3;
use crate::crate_utils::find_deposit_file_name;
use crate::pdb_atom_filters::{ByResidueRange, PdbAtomPredicate};
use crate::pdb_parsing_error::PDBError;

impl Deposit {

    /// Reads PDB-formatted content from a buffer.
    ///
    /// This function allows reading PDB structures from `String`
    /// # Example
    /// ```
    /// use bioshell_pdb::{Deposit};
    /// use std::io::BufReader;
    /// let pdb_txt: &str =
    /// "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
    /// ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C
    /// ATOM     35  CA  TYR A   3      -7.173  -2.314   2.811  1.00  0.08           C
    /// ATOM     56  CA  LYS A   4      -3.922  -3.881   4.044  1.00  0.10           C
    /// ATOM     78  CA  LEU A   5      -0.651  -2.752   2.466  1.00  0.11           C
    /// ATOM     97  CA  ILE A   6       2.338  -5.105   2.255  1.00  0.13           C";
    ///
    /// let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
    /// let seq = deposit.structure().sequence("A");
    /// assert_eq!(seq.to_string(100), "MTYKLI");
    /// ```
    pub fn from_pdb_reader<R: BufRead>(reader: R) -> Result<Deposit, PDBError> {

        let start = Instant::now();

        let mut deposit  = Deposit::new("");

        let mut helices: Vec<PdbHelix> = vec![];
        let mut strands: Vec<PdbSheet> = vec![];
        let mut remarks = PDBRemarks::new();
        let mut seqres: Vec<String> = vec![];
        let mut model_id = 0;
        let mut title : Option<PdbTitle> = None;
        let mut header: Option<PdbHeader> = None;

        let mut pdb_structure = Structure::new("");
        pdb_structure.model_coordinates.push(vec![]);

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
                        if let Some(last_atom) = pdb_structure.atoms.last() {
                            pdb_structure.ter_atoms.insert(last_atom.chain_id.clone(), ResidueId::try_from(last_atom).unwrap());
                        }
                    }
                }
                "HEADER" => {
                    header = PdbHeader::new(&line);
                },
                "EXPDTA" => {
                    deposit.methods = ExperimentalMethod::from_expdata_line(&line);
                },
                "TITLE" => {
                    if title == None {
                        title = Some(PdbTitle::from_pdb_line(&line));
                    } else {
                        let title = title.as_mut().unwrap();
                        title.append_pdb_line(&line);
                    }
                },
                "HELIX" => {
                    let helix = PdbHelix::from_helix_line(&line);
                    helices.push(helix);
                }
                "SHEET" => {
                    let strand = PdbSheet::from_sheet_line(&line);
                    strands.push(strand?);
                }
                "ATOM" | "HETATM" => {
                    if model_id == 0 {
                        let a = PdbAtom::from_atom_line(&line);
                        pdb_structure.model_coordinates[model_id].push(a.pos.clone());
                        pdb_structure.atoms.push(a);
                    } else {
                        let v = Vec3::from_pdb_line(&line);
                        pdb_structure.model_coordinates[model_id].push(v);
                    }
                },
                "ENDMDL" => {
                    model_id += 1;
                },
                "MODEL" => {
                    if pdb_structure.model_coordinates.len() == model_id {
                        pdb_structure.model_coordinates.push(vec![]);
                    }
                },
                "REMARK" => {
                    remarks.add_remark(&line);
                }
                "SEQRES" => {
                    seqres.push(line);
                }
                "CRYST1" => {
                    deposit.unit_cell = Some(UnitCell::from_cryst1_line(&line));
                }
                _ => {},
            };
        }
        debug!("{:} atoms loaded", pdb_structure.atoms.len());

        // ---------- Annotate secondary structure
        for i in 0..helices.len() {
            let from = helices[i].init_res_id();
            let to = helices[i].end_res_id();
            let check = ByResidueRange::new(from,to);
            for a in &mut pdb_structure.atoms {
                if check.check(a) {
                    a.secondary_struct_type
                        = SecondaryStructureTypes::from_pdb_class(helices[i].helix_class as usize)
                }
            }
        }
        for i in 0..strands.len() {
            let from = strands[i].init_res_id();
            let to = strands[i].end_res_id();
            let check = ByResidueRange::new(from,to);
            for a in &mut pdb_structure.atoms {
                if check.check(a) { a.secondary_struct_type = SecondaryStructureTypes::Strand  }
            }
        }

        // ---------- Extract values stored in remarks
        deposit.resolution = remarks.resolution();

        deposit.title = match title {
            None => None,
            Some(title) => {Some(title.to_string())}
        };
        if let Some(header) = header {
            deposit.classification = Some(header.classification);
            deposit.id_code = header.id_code;
            deposit.dep_date = Some(header.dep_date);
        }

        pdb_structure.update();
        deposit.structure = pdb_structure;

        debug!("Structure loaded in: {:?}", start.elapsed());

        Ok(deposit)
    }

    /// Reads a [`Structure`](Structure) from a PDB file
    ///
    pub fn from_pdb_file(file_name: &str) -> Result<Deposit, PDBError> {

        info!("Loading a PDB deposit: {}", file_name);
        let reader = open_file(file_name)?;
        return Self::from_pdb_reader(reader);
    }
}


fn parse_seqres_records(seqres_records: Vec<String>) -> HashMap<String, Sequence> {
    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();

    for record in seqres_records {
        let parts: Vec<&str> = record.split_whitespace().collect();
        if parts.len() < 5 {
            continue; // Skip records that don't have at least one residue entry
        }

        let chain_id = parts[2].to_string();
        let _num_residues: usize = parts[3].parse().unwrap(); // Assuming the number of residues field is correct and parseable

        let residues = &parts[4..]; // The rest are residue names

        // Get the vector for the chain or create it if it doesn't exist
        let chain_sequence = sequences.entry(chain_id).or_insert_with(Vec::new);

        // Add all residues in this record to the vector
        for &residue in residues {
            if let Some(restype) = ResidueTypeManager::get().by_code3(residue) {
                chain_sequence.push(u8::try_from(restype.parent_type.code1()).unwrap());
            } else {
                chain_sequence.push(b'X');
            }
        }
    }
    let mut final_sequences: HashMap<String, Sequence> = HashMap::new();
    for (chain_id, seq) in sequences {
        let header = format!(":{}", chain_id);
        final_sequences.insert(chain_id,  Sequence::from_attrs(header, seq));
    }

    return final_sequences;
}

/// Returns true if a given file is in PDB format.
///
/// This function simply tests whether the first data line of a given file starts with ``HEADER``,
/// ``REMARK``, ``ATOM`` or ``HETATM``.
/// Otherwise, it returns ``false``. When the file can't be open returns I/O error..
pub fn is_pdb_file(file_path: &str) -> io::Result<bool> {
    let reader = open_file(file_path)?;

    let pdb_starts_with = ["HEADER", "ATOM", "HETATM", "REMARK"];
    for line in reader.lines() {
        let line = line?;
        if !line.is_empty() {
            return Ok(pdb_starts_with.iter().any(|s|line.starts_with(s)));
        }
    }

    return Ok(false);
}

static PDB_PREFIXES: [&str; 4] = ["pdb", "PDB", "pdb", ""];
static PDB_SUFFIXES: [&str; 7] = [".ent", ".ent.gz", ".gz", ".pdb", ".PDB", ".pdb.gz", ""];

/// Attempts to find a PDB file in a given directory.
///
/// Looks in the specified path for a file with a given PDB data, identified by
/// a given PDB code. For a given 4-character ID (digit + 3 letters), the method checks
/// the following possibilities:
///
/// - `given_path/1abc`
/// - `given_path/1ABC`
/// - `given_path/1abc.pdb`
/// - `given_path/1ABC.pdb`
/// - `given_path/1ABC.PDB`
/// - `given_path/pdb1abc`
/// - `given_path/PDB1ABC`
/// - `given_path/pdb1abc.ent`
/// - `given_path/PDB1ABC.ent`
/// - `given_path/pdb1abc.ent.gz`
/// - `given_path/PDB1ABC.ent.gz`
/// - `given_path/ab/pdb1abc.ent`
/// - `given_path/ab/pdb1abc.ent.gz`
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
/// use bioshell_pdb::find_pdb_file_name;
/// let result = find_pdb_file_name("2gb1", "./tests/test_files/");
/// assert!(result.is_ok());
/// assert_eq!(result.unwrap(), "./tests/test_files/2gb1.pdb");
/// ```
///
pub fn find_pdb_file_name(pdb_code: &str, pdb_path: &str) -> Result<String, io::Error> {
    find_deposit_file_name(pdb_code, pdb_path, &PDB_PREFIXES, &PDB_SUFFIXES)
}

#[cfg(test)]
mod tests {
    use crate::load_pdb::parse_seqres_records;

    #[test]
    fn test_seqres_loading() {
        let input = vec![
            String::from("SEQRES   1 A   8  ALA VAL CYS LEU MET GLU ARG GLY"),
            String::from("SEQRES   2 A   8  TYR PHE ASN"),
            String::from("SEQRES   1 B   5  LYS THR GLN"),
        ];

        let sequences = parse_seqres_records(input);
        assert_eq!(sequences["A"].to_string(100), "AVCLMERGYFN".to_string());
        assert_eq!(sequences["B"].to_string(100), "KTQ".to_string());
    }
}

