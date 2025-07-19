use std::fmt::{Display, Formatter};

/// The `HEADER` record contains the header data for a PDB entry.
///
/// Typically a [`PdbHeader`](PdbHeader) struct is created by the [`load_pdb_file()`](crate::from_pdb_file()) function
/// when it reads and parses a PDB file.
///
/// See  the [official documentation of the `TITLE` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#HEADER) for details
pub(crate) struct PdbHeader {
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
    pub fn new(header_line: &str) -> Option<Self> {
        Some(Self{
            classification: header_line[10..50].trim().to_string(),
            dep_date: header_line[50..59].trim().to_string(),
            id_code: header_line[59..66].trim().to_string(),
        })
    }
}

impl Display for PdbHeader {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:10}{:40}{:12}{:4}", "HEADER", self.classification, self.dep_date, self.id_code)
    }
}



#[cfg(test)]
mod tests {
    use crate::pdb_header::PdbHeader;

    #[test]
    fn test_pdb_header_new() {
        let header_line = "HEADER    KINASE                                  18-MAY-98   16PK              ".to_string();
        let header = PdbHeader::new(&header_line).unwrap();
        assert_eq!(header.classification, "KINASE".to_string());
        assert_eq!(header.dep_date, "18-MAY-98".to_string());
        assert_eq!(header.id_code, "16PK".to_string());
    }

    #[test]
    fn test_pdb_header_to_string() {
        let header_line = "HEADER    KINASE                                  18-MAY-98   16PK              ".to_string();
        let header = PdbHeader::new(&header_line).unwrap();
        let header_str = header.to_string();
        assert_eq!(header_str, "HEADER    KINASE                                  18-MAY-98   16PK".to_string());
    }
}