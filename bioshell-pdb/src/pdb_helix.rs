use bioshell_cif::{cif_columns_by_name, CifData, parse_item_or_error, CifLoop, CifError, entry_has_value};
use crate::{PDBError, ResidueId, value_or_missing_key_pdb_error};
use crate::PDBError::CifParsingError;
use bioshell_cif::CifError::{ItemParsingError, MissingCifDataKey};

/// Corresponds to a HELIX record.
///
/// HELIX records are used to identify the position of helices in the molecule sequence. The struct
/// provides also a [`ResidueId`](ResidueId) for the initial and the terminal residue of this helix.
///
/// Refer to the [official documentation of the `HELIX` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX)
pub struct PdbHelix {
    /// Serial number of the helix.
    ///
    /// As defined by the PDB standard, it should start at 1  and increases incrementally.
    pub ser_num: i32,
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
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbHelix, ResidueId};
    /// let helix_line = "HELIX    1   A SER A    5  GLN A   11  1                                   7";
    /// let pdb_helix = PdbHelix::from_helix_line(helix_line);
    /// assert_eq!(pdb_helix.init_res_id(), ResidueId::new("A", 5, ' '));
    /// ```
    pub fn init_res_id(&self) -> ResidueId { ResidueId::new(&self.init_chain_id, self.init_seq_num, self.init_i_code)}

    /// Creates [`ResidueId`](ResidueId) struct that identifies the last residue of this helix (inclusive!)
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbHelix, ResidueId};
    /// let helix_line = "HELIX    1   A SER A    5  GLN A   11  1                                   7";
    /// let pdb_helix = PdbHelix::from_helix_line(helix_line);
    /// assert_eq!(pdb_helix.end_res_id(), ResidueId::new("A", 11, ' '));
    /// ```
    pub fn end_res_id(&self) -> ResidueId { ResidueId::new(&self.end_chain_id, self.end_seq_num, self.end_i_code)}

    /// Lists all helices from a given structure
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::PdbHelix;
    /// let cif_data = include_str!("tests/test_files/2gb1.cif");
    /// let reader = BufReader::new(cif_data.as_bytes());
    /// let cif_data = read_cif_buffer(reader).unwrap();
    /// let helices = PdbHelix::from_cif_data(&cif_data[0]).unwrap();
    /// assert_eq!(helices.len(), 1);
    /// let helix = &helices[0];
    /// assert_eq!(helix.length, 15);
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<PdbHelix>, PDBError> {
        fn new_helix(tokens: &[&str]) -> Result<PdbHelix, CifError> {
            Ok(PdbHelix {
                ser_num: 0,
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
                length: 0,
            })
        }
        let mut helices: Vec<PdbHelix> = Vec::new();
        if let Some(helix_loop) = cif_data.first_loop("_struct_conf.id") {
            cif_columns_by_name!(EntityData, "_struct_conf.id",
                "_struct_conf.beg_label_comp_id", "_struct_conf.beg_label_asym_id",
                "_struct_conf.beg_label_seq_id", "_struct_conf.pdbx_beg_PDB_ins_code",
                "_struct_conf.end_label_comp_id", "_struct_conf.end_label_asym_id",
                "_struct_conf.end_label_seq_id", "_struct_conf.pdbx_end_PDB_ins_code",
            );
            let extractor = EntityData::new(helix_loop)?;
            let mut tokens = [""; 9];
            for row in helix_loop.rows() {
                extractor.data_items(&row, &mut tokens);
                let helix = new_helix(&tokens)?;
                helices.push(helix);
            }
            return Ok(helices);
        } else {
            if cif_data.data_items().contains_key("_struct_conf.id") {
                let mut init_i_code = value_or_missing_key_pdb_error!(cif_data, "_struct_conf.pdbx_beg_PDB_ins_code", char);
                if init_i_code == '?' || init_i_code == '.' { init_i_code = ' '; }
                let mut end_i_code = value_or_missing_key_pdb_error!(cif_data, "_struct_conf.pdbx_end_PDB_ins_code", char);
                if end_i_code == '?' || end_i_code == '.' { end_i_code = ' '; }
                helices.push(PdbHelix {
                    ser_num: cif_data.get_item_or_default("_struct_conf.pdbx_PDB_helix_id", 0),
                    helix_id: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.id", String),
                    init_res_name: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.beg_label_comp_id", String),
                    init_chain_id: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.beg_label_asym_id", String),
                    init_seq_num: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.beg_label_seq_id", i32),
                    init_i_code,
                    end_res_name: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.end_label_comp_id", String),
                    end_chain_id: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.end_label_asym_id", String),
                    end_seq_num: value_or_missing_key_pdb_error!(cif_data, "_struct_conf.end_label_seq_id", i32),
                    end_i_code,
                    helix_class: cif_data.get_item_or_default("_struct_conf.pdbx_PDB_helix_class", 0),
                    comment: String::new(),
                    length: cif_data.get_item_or_default("_struct_conf.pdbx_PDB_helix_length", 0),
                });
            }
            Ok(helices)
        }
    }
}