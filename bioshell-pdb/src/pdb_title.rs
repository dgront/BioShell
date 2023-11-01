/// The TITLE record contains a title for a PDB entry.
///
/// Typically a [`PdbTitle`](PdbTitle) struct is created by the [`load_pdb_file()`](crate::load_pdb_file()) function
/// when it reads and parses a PDB file.
///
/// See  the [official documentation of the `TITLE` entry](https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#TITLE) for details
#[derive(PartialEq)]
pub struct PdbTitle { text: String, }

impl PdbTitle {
    /// Create a new record from a given PDB-formatted line.
    ///
    /// The `new()` method accepts only the very first line of the `TITLE` record. Continuation lines
    /// must me appended using the [`append()`](PdbTitle::append()) method.
    pub fn new(line: &str) -> Self {
        Self { text: line[10..].trim().to_string(), }
    }

    /// Appends an additional PDB line to this entry
    ///
    /// According to the PDB file format, the `TITLE` record may be split into multiple lines.
    /// This method allows build the full record by appending lines that follow the first one
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::PdbTitle;
    /// let mut title = PdbTitle::new("TITLE     A NOVEL, HIGHLY STABLE FOLD OF THE IMMUNOGLOBULIN BINDING");
    /// title.append("TITLE    2 DOMAIN OF STREPTOCOCCAL PROTEIN G ");
    /// assert_eq!(title.to_string(),
    ///                    String::from("A NOVEL, HIGHLY STABLE FOLD OF THE IMMUNOGLOBULIN BINDING DOMAIN OF STREPTOCOCCAL PROTEIN G"))
    /// ```
    pub fn append(&mut self, line: &str) {
        self.text.push_str(line[10..].trim_end());
    }

    pub fn to_string(&self) -> String {
        self.text.clone()
    }
}