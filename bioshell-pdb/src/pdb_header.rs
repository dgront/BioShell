use std::fmt::{Display, Formatter};

/// The `HEADER` record contains the header data for a PDB entry.
///
/// Typically a [`PdbHeader`](PdbHeader) struct is created by the [`load_pdb_file()`](crate::load_pdb_file()) function
/// when it reads and parses a PDB file.
///
/// See  the [official documentation of the `TITLE` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#HEADER) for details
///
/// # Example
/// ```
/// use bioshell_pdb::PdbHeader;
/// let header_line = String::from("HEADER    PHOTOSYNTHESIS                          28-MAR-07   2UXK");
/// let h = PdbHeader::new(&header_line);
/// assert_eq!(format!("{}", h), header_line);
/// ```
pub struct PdbHeader {
    /// Classifies the molecule(s)
    ///
    /// This field should contain one of classifications from a curated list available at the [wwPDB website](http://www.wwpdb.org/)
    pub classification: String,
    /// Deposition date
    pub dep_date: String,
    /// Four-character PDB code of this deposit, such as `2GB1` or `4HHB`
    pub id_code: String,
}

impl PdbHeader {
    /// Creates a new [`PdbHeader`](PdbHeader) struct from a PDB-formatted line
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::PdbHeader;
    /// let h = PdbHeader::new("HEADER    PHOTOSYNTHESIS                          28-MAR-07   2UXK");
    /// assert_eq!(h.id_code, String::from("2UXK"));
    /// ```
    pub fn new(header_line: &str) -> Self {
        Self{
            classification: header_line[10..50].trim().to_string(),
            dep_date: header_line[50..59].trim().to_string(),
            id_code: header_line[59..66].trim().to_string(),
        }
    }
}

impl Display for PdbHeader {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:10}{:40}{:12}{:4}", "HEADER", self.classification, self.dep_date, self.id_code)
    }
}