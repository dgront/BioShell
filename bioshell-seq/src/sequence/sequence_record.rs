use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Lines};
use std::path::Path;

/// A general representation of a biological sequence record.
///
/// This struct captures essential metadata and sequence information from
/// flat file-based sources such as GenBank, SwissProt, or similar databases. Each entry provides a sequence
/// (amino acids or nucleotides) and associated metadata, such as the taxid of the source organism,
/// and a unique identifier.
///
/// [`SequenceRecord`] entries can be loaded from files in SwissProt or GenBank formats.
///
/// # Examples
/// ```
/// use bioshell_io::open_file;
/// use bioshell_seq::sequence::SwissProtIterator;
/// use std::io::{BufRead, BufReader, Error};
/// # fn main() -> Result<(), Error> {
/// let reader = open_file("tests/test_files/R4K4X3.spt")?;
/// let mut iterator = SwissProtIterator::new(reader);
/// let record = iterator.next().expect("expected one record");
/// assert_eq!(&record.accession, "R4K4X3");
/// assert_eq!(&record.taxid, &Some(86416));
/// assert_eq!(&record.sequence, "MKGFVDKDTCIGCGLCTSICPEVFIMDDKGKAERSKNEILETLVASAQEAATECPVNAITVE");
/// # Ok(())
/// # }
/// ```
pub struct SequenceRecord {
    /// A unique record identifier, such as a locus name or entry name, depending on the source format.
    pub id: String,

    /// Indicates whether the record has been manually reviewed or curated.
    ///
    /// The interpretation of this flag may vary by source (e.g., `true` for SwissProt entries,
    /// `false` for automated annotations like TrEMBL or GenBank).
    pub is_reviewed: bool,

    /// The accession number that uniquely identifies the record within its database.
    pub accession: String,

    /// The full name of the sequence
    pub full_name: String,

    /// The EC classification code (EC) associated with the sequence, if available.
    pub ec: Option<String>,

    /// The NCBI taxonomy ID (`taxid`) of the source organism, if available.
    pub taxid: Option<u32>,

    /// The scientific name of the species or organism associated with the sequence.
    pub species_name: String,

    /// The raw biological sequence (e.g., nucleotide or amino acid), stored as an uppercase string.
    pub sequence: String,
}

impl SequenceRecord {
    /// Constructs a SwissProt-style FASTA header line (without the leading '>').
    ///
    /// Format: `sp.accession.id Full name OS=Species Name OX=TaxID`
    ///
    /// - `sp` is used if the record is reviewed (`is_reviewed == true`), otherwise `tr`
    /// - `full_name` is always included
    /// - `OS=...` and `OX=...` are included if data is present
    pub fn header(&self) -> String {
        let db_prefix = if self.is_reviewed { "sp" } else { "tr" };
        let mut header = format!("{db_prefix}.{}.{} {}", self.accession, self.id, self.full_name);

        if !self.species_name.is_empty() {
            header.push_str(" OS=");
            header.push_str(&self.species_name);
        }

        if let Some(taxid) = self.taxid {
            header.push_str(&format!(" OX={}", taxid));
        }

        header
    }
}

/// Iterates over [SequenceRecord] structs, parsed from a file ine the **SwissProt** format.
///
/// The Swiss-Prot format provides the amino acid (or nucleic) sequence, gene identifier, accession number,
/// taxid and other information.
/// For details on the  format, see the [UniProt flat file format specification](https://web.expasy.org/docs/userman.html)
///
/// # Example
/// ```
/// use bioshell_io::open_file;
/// use bioshell_seq::sequence::SwissProtIterator;
/// use std::io::{BufRead, BufReader, Error};
/// # fn main() -> Result<(), Error> {
/// let reader = open_file("tests/test_files/A0A068Q5V6.spt")?;
/// let mut iterator = SwissProtIterator::new(reader);
/// let record = iterator.next().expect("expected one record");
/// assert_eq!(&record.accession, "A0A068Q5V6");
/// # Ok(())
/// # }
/// ```
pub struct SwissProtIterator<R: BufRead> {
    lines: Lines<R>,
    current_id: Option<String>,
    current_ac: Option<String>,
    last_de_section: DESectionType,
    current_full: String,
    current_ec: Option<String>,
    current_taxid: Option<u32>,
    current_species: String,
    current_sequence: String,
    current_is_reviewed: bool,
    expected_sequence_length: usize,
}

impl<R: BufRead> SwissProtIterator<R> {
    pub fn new(reader: R) -> Self {
        SwissProtIterator {
            lines: reader.lines(),
            current_id: None,
            current_ac: None,
            last_de_section: DESectionType::None,
            current_full: String::new(),
            current_ec: None,
            current_taxid: None,
            current_species: String::new(),
            current_sequence: String::new(),
            current_is_reviewed: false,
            expected_sequence_length: 0,
        }
    }

    fn reset(&mut self) {
        self.current_id = None;
        self.current_ac = None;
        self.last_de_section = DESectionType::None;
        self.current_full.clear();
        self.current_ec = None;
        self.current_taxid = None;
        self.current_species.clear();
        self.current_sequence.clear();
        self.current_is_reviewed = false;
        self.expected_sequence_length = 0;
    }

    fn parse_id(&mut self, line: &str) {
        let tokens: Vec<&str> = line[5..].split_whitespace().collect();

        if tokens.len() >= 3 {
            self.current_id = Some(tokens[0].to_string());
            self.current_is_reviewed = matches!(tokens[1].trim_end_matches(';'), "Reviewed");

            // Parse length (e.g., "104" from "104 AA.", the "AA" string is not required)
            if let Ok(length) = tokens[2].parse::<usize>() {
                self.expected_sequence_length = length;
            }
        }
    }

    fn parse_ac(&mut self, line: &str) {
        if self.current_ac.is_none() {
            if let Some(first) = line[5..].split(';').map(str::trim).find(|s| !s.is_empty()) {
                self.current_ac = Some(first.to_string());
            } else {
                panic!("An entry must have only one unique accession !")
            }
        }
    }

    /// Parses a `DE` line and extracts only the `Full=` and `EC=` fields from the `RecName` section.
    pub fn parse_de(&mut self, line: &str) {
        let trimmed = line.trim_start();

        if trimmed.starts_with("DE   RecName:") {
            self.last_de_section = DESectionType::RecName;
            self.parse_de_fields(&trimmed["DE   RecName:".len()..].trim_start());
        } else if trimmed.starts_with("DE   AltName:") {
            self.last_de_section = DESectionType::AltName;
        } else if trimmed.starts_with("DE   SubName:") {
            self.last_de_section = DESectionType::SubName;
        } else if trimmed.starts_with("DE") && self.last_de_section == DESectionType::RecName {
            self.parse_de_fields(&trimmed["DE".len()..].trim_start());
        } else {
            self.last_de_section = DESectionType::None;
        }
    }

    /// Parses `Full=` or `EC=` fields from a line assumed to belong to the `RecName` section.
    fn parse_de_fields(&mut self, content: &str) {
        if let Some(start) = content.find("Full=") {
            let after_eq = &content[start + "Full=".len()..];
            let end = after_eq.find(|c| c == ';' || c == '{').unwrap_or(after_eq.len());
            let value = after_eq[..end].trim();
            if !value.is_empty() {
                self.current_full = value.to_string();
            }
        }

        if let Some(start) = content.find("EC=") {
            let after_eq = &content[start + "EC=".len()..];
            let end = after_eq.find(|c| c == ';' || c == '{').unwrap_or(after_eq.len());
            let value = after_eq[..end].trim();
            if !value.is_empty() {
                self.current_ec = Some(value.to_string());
            }
        }
    }

    /// Parse the OS (Organism Species) line
    fn parse_os(&mut self, line: &str) {
        // Append line content without prefix, preserving space between lines
        let content = line[5..].trim();
        if !self.current_species.is_empty() {
            self.current_species.push(' ');
        }
        self.current_species.push_str(content);
    }

    /// Parse the OX (Organism taxonomy cross-reference) line
    ///
    /// # Examples of valid input lines:
    /// - `OX   NCBI_TaxID=1229662;`
    /// - `OX NCBI_TaxID=4058 {ECO:0000312|EMBL:AEB69788.1};`
    /// - `OX   NCBI_TaxID=86416 {ECO:0000313|EMBL:AGK95569.1, ECO:0000313|Proteomes:UP000013523};`
    fn parse_ox(&mut self, line: &str) {
        for token in line.split_whitespace() {
            if let Some(taxid_str) = token.strip_prefix("NCBI_TaxID=") {
                let end = taxid_str
                    .find(|c: char| !c.is_ascii_digit())
                    .unwrap_or_else(|| taxid_str.len());

                let numeric_part = &taxid_str[..end];

                if let Ok(parsed) = numeric_part.parse::<u32>() {
                    self.current_taxid = Some(parsed);
                }
                break;
            }
        }
    }

    fn parse_sq_line(&mut self, line: &str) {
        self.current_sequence.push_str(&line.replace(' ', "").trim());
    }
}

macro_rules! finalize_record {
    ($self:ident) => {
        if let (Some(id), Some(ac)) = (&$self.current_id, &$self.current_ac) {
            let sequence = std::mem::take(&mut $self.current_sequence);
            if let expected_len = $self.expected_sequence_length {
                if sequence.len() != expected_len {
                    eprintln!(
                        "Warning: sequence length mismatch for entry {}: expected {}, got {}",
                        id, expected_len, sequence.len()
                    );
                }
            }
            let record = SequenceRecord {
                id: id.clone(),
                is_reviewed: $self.current_is_reviewed,
                accession: ac.clone(),
                full_name: std::mem::take(&mut $self.current_full),
                ec: std::mem::take(&mut $self.current_ec),
                taxid: $self.current_taxid,
                species_name: std::mem::take(&mut $self.current_species),
                sequence,
            };
            $self.reset();
            return Some(record);
        }
    };
}

impl<R: BufRead> Iterator for SwissProtIterator<R> {
    type Item = SequenceRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut in_sequence = false;

        while let Some(Ok(line)) = self.lines.next() {
            if line.starts_with("//") {
                finalize_record!(self);
            } else if line.starts_with("ID") {
                self.parse_id(&line);
            } else if line.starts_with("AC") {
                self.parse_ac(&line);
            } else if line.starts_with("DE") {
                self.parse_de(&line);
            } else if line.starts_with("OX") {
                self.parse_ox(&line);
            } else if line.starts_with("OS") {
                self.parse_os(&line);
            } else if line.starts_with("SQ") {
                in_sequence = true;
            } else if in_sequence {
                self.parse_sq_line(&line);
            }
        }

        // Handle case where file ends without a final `//`
        finalize_record!(self);

        None
    }
}

#[derive(Debug, PartialEq)]
enum DESectionType {
    RecName,
    AltName,
    SubName,
    None,
}

/// Iterates over [SequenceRecord] structs, parsed from all the SwissProt files found in a given folder.
///
/// [SwissProtFolderIterator] allows convenient processing of all [SequenceRecord] structs found in
/// all files in a given folder.
/// The Swiss-Prot format provides the amino acid (or nucleic) sequence, gene identifier, accession number,
/// taxid and other information.
/// For details on the  format, see the [UniProt flat file format specification](https://web.expasy.org/docs/userman.html)
///
/// # Example
/// ```
/// use bioshell_io::open_file;
/// use bioshell_seq::sequence::SwissProtFolderIterator;
/// use std::io::{Error};
/// # fn main() -> Result<(), Error> {
/// let mut iterator = SwissProtFolderIterator::from_folder("tests/test_files/", "spt")?;
/// let mut accession_found: Vec<String> = vec![];
/// for seq_record in iterator {
///     accession_found.push(seq_record.accession)
/// }
/// # assert_eq!(accession_found.len(), 2);
/// # Ok(())
/// # }
/// ```
pub struct SwissProtFolderIterator {
    inner: Box<dyn Iterator<Item = SequenceRecord>>,
}

impl SwissProtFolderIterator {
    pub fn from_folder<P: AsRef<Path>>(folder: P, extension: &str) -> std::io::Result<Self> {
        let mut iterators: Vec<Box<dyn Iterator<Item = SequenceRecord>>> = Vec::new();

        for entry in fs::read_dir(folder)? {
            let entry = entry?;
            let path = entry.path();

            if path.extension().map_or(false, |ext| ext == extension) {
                let file = File::open(&path)?;
                let reader = BufReader::new(file);
                let sp_iter = SwissProtIterator::new(reader);
                iterators.push(Box::new(sp_iter));
            }
        }

        // Flatten into one iterator
        let unified = iterators.into_iter().flatten();
        Ok(Self {
            inner: Box::new(unified),
        })
    }
}

impl Iterator for SwissProtFolderIterator {
    type Item = SequenceRecord;

    fn next(&mut self) -> Option<Self::Item> { self.inner.next() }
}


/// Iterates over [SequenceRecord] structs, parsed from a file ine the **NCBI / GeneBank** format.
///
/// The NCBI format provides the amino acid (or nucleic) sequence, gene identifier, accession number,
/// taxid and other information.
    /// For details on the  format, see the [Gene Bank flat file format specification](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
///
/// # Example
/// ```
/// use bioshell_io::open_file;
/// use bioshell_seq::sequence::NCBIIterator;
/// use std::io::{BufRead, BufReader, Error};
/// # fn main() -> Result<(), Error> {
/// let reader = open_file("tests/test_files/WP_015613896.gp")?;
/// let mut iterator = NCBIIterator::new(reader);
/// let record = iterator.next().expect("expected one record");
/// assert_eq!(&record.accession, "WP_015613896");
/// # Ok(())
/// # }
/// ```
pub struct NCBIIterator<R: BufRead> {
    lines: Lines<R>,
}

impl<R: BufRead> NCBIIterator<R> {
    pub fn new(reader: R) -> Self {
        Self {
            lines: reader.lines(),
        }
    }
}

impl<R: BufRead> Iterator for NCBIIterator<R> {
    type Item = SequenceRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut accession = String::new();
        let mut id = String::new();
        let mut species_name = String::new();
        let mut taxid = None;
        let mut sequence = String::new();
        let mut in_origin = false;

        while let Some(Ok(line)) = self.lines.next() {
            let trimmed = line.trim_start();

            if trimmed.starts_with("LOCUS") {
                id = trimmed.split_whitespace().nth(1).unwrap_or("").to_string();
            } else if trimmed.starts_with("ACCESSION") {
                accession = trimmed.split_whitespace().nth(1).unwrap_or("").to_string();
            } else if trimmed.starts_with("SOURCE") {
                species_name = trimmed["SOURCE".len()..].trim().to_string();
            } else if trimmed.starts_with("/db_xref=\"taxon:") {
                // e.g. /db_xref="taxon:9606"
                if let Some(taxid_str) = trimmed.split("taxon:").nth(1) {
                    if let Some(end) = taxid_str.find('"') {
                        taxid = taxid_str[..end].parse().ok();
                    }
                }
            } else if trimmed.starts_with("ORIGIN") {
                in_origin = true;
            } else if trimmed.starts_with("//") {
                // End of record
                return Some(SequenceRecord {
                    id,
                    is_reviewed: false,
                    accession,
                    full_name: "".to_string(),
                    ec: None,
                    taxid,
                    species_name,
                    sequence,
                });
            } else if in_origin {
                // Sequence lines: "     1 atgccc..."
                let seq_part: String = trimmed
                    .split_whitespace()
                    .skip(1)
                    .collect::<Vec<&str>>()
                    .join("")
                    .to_ascii_uppercase(); // ← convert to uppercase here
                sequence.push_str(&seq_part);
            }
        }

        // End of file: if at least a partial record was read
        if !accession.is_empty() {
            Some(SequenceRecord {
                id,
                is_reviewed: false,
                accession,
                full_name: "".to_string(),
                ec: None,
                taxid,
                species_name,
                sequence,
            })
        } else {
            None
        }
    }
}
