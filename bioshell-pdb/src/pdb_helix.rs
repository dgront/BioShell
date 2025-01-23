use bioshell_cif::{CifData, parse_item_or_error, CifError, entry_has_value, CifTable};
use crate::{PDBError, ResidueId};
use bioshell_cif::CifError::{ItemParsingError};

/// Corresponds to a HELIX record.
///
/// HELIX records are used to identify the position of helices in the molecule sequence. The struct
/// provides also a [`ResidueId`](ResidueId) for the initial and the terminal residue of this helix.
///
/// Refer to the [official documentation of the `HELIX` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX)
#[allow(dead_code)]
pub(crate) struct PdbHelix {

    /// Unique identifier for the helix.
    ///
    /// In a PDB file it is an integer number, but in a CIF file it is a string.
    pub ser_num: String,
    /// Helix  identifier.
    ///
    /// In addition to a serial number, each helix is given an alphanumeric character helix identifier.
    pub helix_id: String,
    /// Name of the initial residue in the 3-letter code
    pub init_res_name: String,
    /// Chain identifier for the initial residue
    pub init_chain_id: String,
    /// Sequence number of the initial residue
    pub init_seq_num: i32,
    /// Insertion code of the initial residue
    pub init_i_code: char,
    /// Name of the terminal residue of the helix
    pub end_res_name: String,
    /// Chain identifier for the terminal residue
    pub end_chain_id: String,
    /// Sequence number of the terminal residue
    pub end_seq_num: i32,
    /// Insertion code of the terminal residue
    pub end_i_code: char,
    /// Helix class: an integer from 1 to 12 (both inclusive)
    ///
    /// For the list of codes see the [official documentation of the `HELIX` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX)
    pub helix_class: u8,
    /// Comment about this helix
    pub comment: String,
    /// Length of this helix
    pub length: i32,
}

impl PdbHelix {

    /// Creates a [`PdbHelix`](PdbHelix) struct from a relevant PDB-formatted line.
    pub fn from_helix_line(line: &str) -> PdbHelix {

        let ser_num = line[7..10].trim().parse().unwrap();
        let helix_id = line[11..14].trim().to_string();
        let init_res_name = line[15..18].trim().to_string();
        let init_chain_id = line[19..20].to_string();
        let init_seq_num = line[21..25].trim().parse().unwrap();
        let init_i_code = line.chars().nth(26).unwrap();
        let end_res_name = line[27..30].trim().to_string();
        let end_chain_id = line[31..32].to_string();
        let end_seq_num = line[33..37].trim().parse().unwrap();
        let end_i_code = line.chars().nth(37).unwrap();
        let helix_class = line[38..40].trim().parse().unwrap();
        let comment = line[40..70].to_string();
        let length = line[71..76].trim().parse().unwrap();
        PdbHelix { ser_num, helix_id, init_res_name, init_chain_id, init_seq_num, init_i_code,
            end_res_name, end_chain_id, end_seq_num, end_i_code, helix_class, comment, length, }
    }

    /// Creates [`ResidueId`](ResidueId) struct that identifies the first residue of this helix
    pub fn init_res_id(&self) -> ResidueId { ResidueId::new(&self.init_chain_id, self.init_seq_num, self.init_i_code)}

    /// Creates [`ResidueId`](ResidueId) struct that identifies the last residue of this helix (inclusive!)
    pub fn end_res_id(&self) -> ResidueId { ResidueId::new(&self.end_chain_id, self.end_seq_num, self.end_i_code)}

    /// Lists all helices from a given structure.
    ///
    /// Note, that the returned vector may be empty if the structure does not contain any helix.
    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<PdbHelix>, PDBError> {
        fn new_helix(tokens: &[&str]) -> Result<PdbHelix, CifError> {
            Ok(PdbHelix {
                ser_num: "1".to_string(),
                helix_id: tokens[0].to_string(),
                init_res_name: tokens[1].to_string(),
                init_chain_id: tokens[2].to_string(),
                init_seq_num: parse_item_or_error!(tokens[3], i32),
                init_i_code: if !entry_has_value(tokens[4]) { ' ' } else { tokens[4].chars().nth(0).unwrap() },
                end_res_name: tokens[5].to_string(),
                end_chain_id: tokens[6].to_string(),
                end_seq_num: parse_item_or_error!(tokens[7], i32),
                end_i_code: if !entry_has_value(tokens[4]) { ' ' } else { tokens[8].chars().nth(0).unwrap() },
                helix_class: 1,
                comment: String::new(),
                length: parse_item_or_error!(tokens[9], i32),
            })
        }
        let mut helices: Vec<PdbHelix> = Vec::new();
        if let Ok(helix_table) = CifTable::new(cif_data, "_struct_conf.",
["id",  "beg_auth_comp_id", "beg_auth_asym_id", "beg_auth_seq_id",
                        "pdbx_beg_PDB_ins_code", "end_auth_comp_id", "end_auth_asym_id",
                        "end_auth_seq_id", "pdbx_end_PDB_ins_code", "pdbx_PDB_helix_length", ]) {
            for row in helix_table.iter() {
                let helix = new_helix(&row)?;
                helices.push(helix);
            }
        }

        return Ok(helices);
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use crate::pdb_helix::PdbHelix;
    use crate::{PDBError, ResidueId};

    #[test]
    fn test_helix_from_line() {
        let helix_line = "HELIX    1   A SER A    5  GLN A   11  1                                   7";
        let pdb_helix = PdbHelix::from_helix_line(helix_line);
        assert_eq!(pdb_helix.init_res_id(), ResidueId::new("A", 5, ' '));
        assert_eq!(pdb_helix.end_res_id(), ResidueId::new("A", 11, ' '));
    }


    #[test]
    fn helix_from_2gb1() -> Result<(), PDBError> {
        let cif_data = include_str!("../tests/test_files/2gb1.cif");
        let reader = BufReader::new(cif_data.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(helices.len(), 1);
        let helix = &helices[0];
        assert_eq!(helix.length, 15);

        Ok(())
   }
    #[test]
    fn helices_from_cif() {
        #[allow(non_upper_case_globals)]
        const cif_2gb1: &str = include_str!("../tests/test_files/2gb1.cif");

        let reader = BufReader::new(cif_2gb1.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(helices.len(), 1);
        let helix = &helices[0];
        assert_eq!(helix.length, 15);
        assert_eq!(helix.ser_num, "1");
        assert_eq!(helix.init_chain_id, "A");
        assert_eq!(helix.end_chain_id, "A");
        assert_eq!(helix.init_seq_num, 22);
        assert_eq!(helix.end_seq_num, 36);
        assert_eq!(helix.init_i_code, ' ');
        assert_eq!(helix.end_i_code, ' ');

        #[allow(non_upper_case_globals)]
        const cif_2fdo: &str = include_str!("../tests/test_files/2fdo.cif");

        let reader = BufReader::new(cif_2fdo.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(helices.len(), 8);
    }
}