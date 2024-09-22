use bioshell_cif::{cif_columns_by_name, CifData, CifError, CifLoop, entry_has_value, parse_item_or_error};
use bioshell_cif::CifError::{MissingCifLoopKey, ItemParsingError};
use crate::{PDBError, ResidueId};
use crate::PDBError::CifParsingError;

/// Corresponds to a SHEET record.
///
/// SHEET records are used to identify the position of beta-strands in the molecule sequence. The struct
/// provides also a [`ResidueId`](ResidueId) for the initial and the terminal residue of this strand.
///
/// Refer to the [official documentation of the `SHEET` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET)
pub struct PdbSheet {
    /// strand  number
    ///
    /// starts at 1 for each strand within a sheet
    pub strand: i32,
    /// sheet  identifier
    pub sheet_id: String,
    /// number  of strands in sheet
    pub num_strands: i32,
    /// name of the initial residue for this strand
    pub init_res_name: String,
    /// chain identifier of the initial residue for this strand
    pub init_chain_id: String,
    /// sequence number of the initial residue  for this strand
    pub init_seq_num: i32,
    /// insertion code of the initial residue  for this strand
    pub init_i_code: char,
    /// name of the terminal residue for this strand
    pub end_res_name: String,
    /// chain identifier of the terminal residue for this strand
    pub end_chain_id: String,
    /// sequence number of the terminal residue  for this strand
    pub end_seq_num: i32,
    /// insertion code of the terminal residue  for this strand
    pub end_i_code: char,
    /// sense of strand with respect to previous strand in the sheet.
    ///
    /// 0 if first strand, 1 if  parallel,and -1 if anti-parallel
    pub sense: i8,
}

impl PdbSheet {

    /// Creates a [`PdbSheet`](PdbSheet) struct from a relevant PDB-formatted line.
    pub fn from_sheet_line(line: &str) -> PdbSheet {

        let strand = line[7..10].trim().parse().unwrap();
        let sheet_id = line[11..14].trim().to_string();
        let num_strands = line[14..16].trim().parse().unwrap();
        let init_res_name = line[17..20].trim().to_string();
        let init_chain_id = line[21..22].to_string();
        let init_seq_num = line[22..26].trim().parse().unwrap();
        let init_i_code = line.chars().nth(26).unwrap();

        let end_res_name = line[28..31].trim().to_string();
        let end_chain_id = line[32..33].to_string();
        let end_seq_num = line[33..37].trim().parse().unwrap();
        let end_i_code = line.chars().nth(37).unwrap();
        let sense = line[38..40].trim().parse().unwrap();
        PdbSheet { strand, sheet_id, num_strands, init_res_name, init_chain_id, init_seq_num, init_i_code,
            end_res_name, end_chain_id, end_seq_num, end_i_code, sense }
    }

    /// Creates [`ResidueId`](ResidueId) struct that identifies the first residue of this helix
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbSheet, ResidueId};
    /// let sheet_line = "SHEET    1   A 5 THR A 107  ARG A 110  0";
    /// let pdb_strand = PdbSheet::from_sheet_line(sheet_line);
    /// assert_eq!(pdb_strand.init_res_id(), ResidueId::new("A", 107, ' '));
    /// ```
    pub fn init_res_id(&self) -> ResidueId { ResidueId::new(&self.init_chain_id, self.init_seq_num, self.init_i_code)}

    /// Creates [`ResidueId`](ResidueId) struct that identifies the last residue of this helix (inclusive!)
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbSheet, ResidueId};
    /// let helix_line = "SHEET    1   A 5 THR A 107  ARG A 110  0";
    /// let pdb_helix = PdbSheet::from_sheet_line(helix_line);
    /// assert_eq!(pdb_helix.end_res_id(), ResidueId::new("A", 110, ' '));
    /// ```
    pub fn end_res_id(&self) -> ResidueId { ResidueId::new(&self.end_chain_id, self.end_seq_num, self.end_i_code)}

    /// Lists all beta-sheets found in a given structure
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::PdbSheet;
    /// let cif_data = include_str!("../tests/test_files/2gb1.cif");
    /// let reader = BufReader::new(cif_data.as_bytes());
    /// let cif_data = read_cif_buffer(reader).unwrap();
    /// let sheets = PdbSheet::from_cif_data(&cif_data[0]).unwrap();
    /// assert_eq!(sheets.len(), 4);
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<PdbSheet>, PDBError> {
        fn new_sheet(tokens: &[&str]) -> Result<PdbSheet, CifError> {
            Ok(PdbSheet {
                strand: parse_item_or_error!(tokens[1], i32),
                sheet_id: tokens[0].to_string(),
                num_strands: 0,
                init_res_name: tokens[2].to_string(),
                init_chain_id: tokens[3].to_string(),
                init_seq_num: parse_item_or_error!(tokens[4], i32),
                init_i_code: if !entry_has_value(tokens[5]) { ' ' } else { tokens[5].chars().nth(0).unwrap() },
                end_res_name: tokens[6].to_string(),
                end_chain_id: tokens[7].to_string(),
                end_seq_num: parse_item_or_error!(tokens[8], i32),
                end_i_code: if !entry_has_value(tokens[9]) { ' ' } else { tokens[9].chars().nth(0).unwrap() },
                sense: 0,
            })
        }
        let mut strands: Vec<PdbSheet> = Vec::new();
        if let Some(strands_loop) = cif_data.first_loop("_struct_sheet_range.id") {
            cif_columns_by_name!(EntityData, "_struct_sheet_range.sheet_id","_struct_sheet_range.id",
                "_struct_sheet_range.beg_label_comp_id","_struct_sheet_range.beg_label_asym_id",
                "_struct_sheet_range.beg_label_seq_id","_struct_sheet_range.pdbx_beg_PDB_ins_code",
                "_struct_sheet_range.end_label_comp_id","_struct_sheet_range.end_label_asym_id",
                "_struct_sheet_range.end_label_seq_id","_struct_sheet_range.pdbx_end_PDB_ins_code",
            );
            let extractor = EntityData::new(strands_loop)?;
            let mut tokens = [""; 10];
            for row in strands_loop.rows() {
                extractor.data_items(&row, &mut tokens);
                let helix = new_sheet(&tokens)?;
                strands.push(helix);
            }
            return Ok(strands);
        } else { Err(CifParsingError{ 0: MissingCifLoopKey {item_key: "_struct_sheet_range.id".to_string()} } ) }
    }
}