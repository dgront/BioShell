use crate::ResidueId;

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
}