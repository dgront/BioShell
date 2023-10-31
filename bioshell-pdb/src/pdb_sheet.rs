use crate::ResidueId;

/// Corresponds to a SHEET record.
///
/// SHEET records are used to identify the position of beta-strands in the molecule sequence. The struct
/// provides also a [`ResidueId`](ResidueId) for the initial and the terminal residue of this strand.
///
/// Refer to the [official documentation of the `SHEET` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#SHEET)
pub struct PdbSheet {
    pub strand: i32,
    pub sheet_id: String,
    pub num_strands: i32,
    pub init_res_name: String,
    pub init_chain_id: String,
    pub init_seq_num: i32,
    pub init_i_code: char,
    pub end_res_name: String,
    pub end_chain_id: String,
    pub end_seq_num: i32,
    pub end_i_code: char,
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
    /// let helix_line = "SHEET    1   A 5 THR A 107  ARG A 110  0";
    /// let pdb_strand = PdbSheet::from_sheet_line(helix_line);
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

}