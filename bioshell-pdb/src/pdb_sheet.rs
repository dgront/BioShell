use bioshell_cif::{CifData, CifError, CifTable, entry_has_value, parse_item_or_error};
use bioshell_cif::CifError::{ItemParsingError};
use crate::{PDBError, ResidueId};

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
}

impl PdbSheet {

    /// Creates a [`PdbSheet`](PdbSheet) struct from a relevant PDB-formatted line.
    pub fn from_sheet_line(line: &str) -> Result<PdbSheet, PDBError> {

        let strand = line[7..10].trim().parse().map_err(|_| PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;
        let sheet_id = line[11..14].trim().to_string();
        let num_strands = line[14..16].trim().parse().map_err(|_| PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;
        let init_res_name = line[17..20].trim().to_string();
        let init_chain_id = line[21..22].to_string();
        let init_seq_num = line[22..26].trim().parse().map_err(|_| PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;
        let init_i_code = line.chars().nth(26).ok_or(PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;

        let end_res_name = line[28..31].trim().to_string();
        let end_chain_id = line[32..33].to_string();
        let end_seq_num = line[33..37].trim().parse().map_err(|_| PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;
        let end_i_code = line.chars().nth(37).ok_or(PDBError::InvalidPdbLineFormat { broken_pdb_line: line.to_string() })?;

        Ok(PdbSheet { strand, sheet_id, num_strands, init_res_name, init_chain_id, init_seq_num, init_i_code,
            end_res_name, end_chain_id, end_seq_num, end_i_code })
    }

    /// Creates [`ResidueId`](ResidueId) struct that identifies the first residue of this helix
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbSheet, ResidueId, PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// let sheet_line = "SHEET    1   A 5 THR A 107  ARG A 110  0";
    /// let pdb_strand = PdbSheet::from_sheet_line(sheet_line)?;
    /// assert_eq!(pdb_strand.init_res_id(), ResidueId::new("A", 107, ' '));
    /// # Ok(())
    /// # }
    /// ```
    pub fn init_res_id(&self) -> ResidueId { ResidueId::new(&self.init_chain_id, self.init_seq_num, self.init_i_code)}

    /// Creates [`ResidueId`](ResidueId) struct that identifies the last residue of this helix (inclusive!)
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PdbSheet, ResidueId, PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// let strand_line = "SHEET    1   A 5 THR A 107  ARG A 110  0";
    /// let pdb_strand = PdbSheet::from_sheet_line(strand_line)?;
    /// assert_eq!(pdb_strand.end_res_id(), ResidueId::new("A", 110, ' '));
    /// # Ok(())
    /// # }
    /// ```
    pub fn end_res_id(&self) -> ResidueId { ResidueId::new(&self.end_chain_id, self.end_seq_num, self.end_i_code)}

    /// Lists all beta-sheets found in a given structure
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::{PdbSheet,PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/2gb1.cif");
    /// let reader = BufReader::new(cif_data.as_bytes());
    /// let cif_data = read_cif_buffer(reader).unwrap();
    /// let sheets = PdbSheet::from_cif_data(&cif_data[0])?;
    /// assert_eq!(sheets.len(), 4);
    /// # Ok(())
    /// # }
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
            })
        }
        let mut strands: Vec<PdbSheet> = Vec::new();
        let strands_table = CifTable::new(cif_data, "_struct_sheet_range",
            [".sheet_id", ".id", ".beg_label_comp_id",".beg_label_asym_id", ".beg_label_seq_id",
            ".pdbx_beg_PDB_ins_code", ".end_label_comp_id",".end_label_asym_id",
            ".end_label_seq_id",".pdbx_end_PDB_ins_code",]
        )?;
        for tokens in strands_table.iter() {
            let strand = new_sheet(&tokens)?;
            strands.push(strand);
        }
        return Ok(strands);
    }
}