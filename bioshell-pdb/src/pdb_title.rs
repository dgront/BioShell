/// The `TITLE` record contains a title for a PDB entry.
///
/// Typically a [`PdbTitle`](PdbTitle) struct is created by the [`load_pdb_file()`](crate::load_pdb_file()) function
/// when it reads and parses a PDB file.
///
/// See  the [official documentation of the `TITLE` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#TITLE) for details
#[derive(PartialEq)]
pub struct PdbTitle { pub text: String, }

impl PdbTitle {

    /// Create a new record from a given PDB-formatted line.
    ///
    /// The `new()` method accepts only the very first line of the `TITLE` record. Continuation lines
    /// must me appended using the [`append_pdb_line()`](PdbTitle::append_pdb_line()) method.
    pub fn from_pdb_line(line: &str) -> Self {
        Self { text: line[10..].trim().to_string(), }
    }

    /// Appends an additional PDB line to this entry
    ///
    /// According to the PDB file format, the `TITLE` record may be split into multiple lines.
    /// This method allows build the full record by appending lines that follow the first one
    pub fn append_pdb_line(&mut self, line: &str) {
        self.text.push_str(line[10..].trim_end());
    }

    pub fn to_string(&self) -> String {
        self.text.clone()
    }
}

#[cfg(test)]
mod tests {
    use crate::pdb_title::PdbTitle;

    #[test]
    fn test_pdb_title() {
        let line = "TITLE     PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE           ";
        let expected_text = "PHOSPHOGLYCERATE KINASE FROM TRYPANOSOMA BRUCEI BISUBSTRATE";
        let pdb_title = PdbTitle::from_pdb_line(line);

        assert_eq!(pdb_title.to_string(), expected_text);
    }
}