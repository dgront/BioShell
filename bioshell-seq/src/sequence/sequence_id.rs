use std::cmp::Ordering;
use std::fmt;
use std::ops::{Deref, DerefMut};
use regex::Regex;

use bioshell_io::sanitize_filename;

/// Represents biological sequence identifiers from major databases.
///
/// This enum captures various standard sequence ID types used in molecular biology and bioinformatics.
/// Each variant holds the original matched identifier as a `String`.
///
/// ```
/// use bioshell_seq::sequence::SeqId;
/// let id = SeqId::RefSeq("XP_123456.1".to_string());
/// match id {
///     SeqId::RefSeq(accession) => println!("RefSeq ID: {}", accession),
///     _ => println!("Other ID type"),
/// }
/// ```
///
/// The [`parse_sequence_id`]() function can be used to parse a sequence description string into a list of `SeqId` variants.
/// For convenience,  [`parse_sequence_id`]() returns [`SeqIdList`].
/// ```
/// use bioshell_seq::sequence::{parse_sequence_id, SeqId};
/// let ids = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9 2' cyclic ADP-D-ribose synthase [taxid=1310613]");
/// assert_eq!(ids.len(), 3);
/// assert!(matches!(ids[0], SeqId::UniProtID(_)));
/// assert!(matches!(ids[1], SeqId::TrEmbl(_)));
/// assert!(matches!(ids[2], SeqId::TaxId(_)));
/// assert_eq!(&ids.to_string(), "UniProtID|ABTIR_ACIB9|tr|A0A009IHW8|[taxid=1310613]");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SeqId {
    /// PDB (Protein Data Bank) 4-character structure ID, optionally with chain (e.g., "1HHP", "1HHP:A").
    PDB(String),

    /// SwissProt or UniProtKB/TrEMBL accession number (e.g., "P12345", "Q9NQX5-2").
    SwissProt(String),

    /// UniProt entry name (e.g., "P53_HUMAN", "RL21_YEAST").
    UniProtID(String),

    /// UniRef (UniProt Reference Clusters), e.g., "UniRef100_P12345".
    UniRef(String),

    /// RefSeq protein or transcript accession from NCBI (e.g., "XP_123456.1", "NM_001256789.2").
    RefSeq(String),

    /// GenBank or EMBL accession (e.g., "AB123456.1", "U49845").
    GenBank(String),

    /// Ensembl gene, transcript, or protein ID (e.g., "ENSG00000139618", "ENST00000331789").
    Ensembl(String),

    /// TrEmbl section of UniProtKB (e.g., "tr|A0A009IHW8|").
    TrEmbl(String),

    /// NCBI GI number ("gi|12345678" or "GI:12345678").
    NCBIGI(String),

    /// NCBI Taxonomy ID (e.g., "[taxid=9606]", "[taxid=10090]").
    TaxId(String),

    /// If nothing has been found, use the first word of the description
    Default(String)
}

impl PartialOrd for SeqId {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.order_priority().cmp(&other.order_priority()))
    }
}

impl Ord for SeqId {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl SeqId {
    fn order_priority(&self) -> usize {
        match self {
            SeqId::PDB(_) => 0,
            SeqId::SwissProt(_) => 1,
            SeqId::UniProtID(_) => 2,
            SeqId::UniRef(_) => 3,
            SeqId::RefSeq(_) => 4,
            SeqId::GenBank(_) => 5,
            SeqId::Ensembl(_) => 6,
            SeqId::TrEmbl(_) => 7,
            SeqId::NCBIGI(_) => 8,
            SeqId::Default(_) => 10,
            SeqId::TaxId(_) => 11,
        }
    }
}

impl fmt::Display for SeqId {
    /// Formats the sequence ID for display.
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::sequence::SeqId;
    /// let seq_id = SeqId::RefSeq("XP_001234567.1".to_string());
    /// let header = seq_id.to_string();
    /// assert_eq!(header, "RefSeq|XP_001234567.1");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SeqId::PDB(s) => write!(f, "pdb|{}", s),
            SeqId::SwissProt(s) => write!(f, "sp|{}", s),
            SeqId::UniProtID(s) => write!(f, "UniProtID|{}", s),
            SeqId::UniRef(s) => write!(f, "UniRef|{}", s),
            SeqId::RefSeq(s) => write!(f, "RefSeq|{}", s),
            SeqId::GenBank(s) => write!(f, "gb|{}", s),
            SeqId::Ensembl(s) => write!(f, "Ensembl|{}", s),
            SeqId::NCBIGI(s) => write!(f, "gi|{}", s),
            SeqId::TrEmbl(s) => write!(f, "tr|{}", s),
            SeqId::TaxId(s) => write!(f, "[taxid={}]", s),
            SeqId::Default(s) => write!(f, "{}", s),
        }
    }
}

/// Extract sequence identifiers from a free-text description string.
///
/// This function scans the input for standard database identifiers such as those from
/// PDB, SwissProt/TrEMBL, UniRef, RefSeq, GenBank/EMBL, Ensembl, and NCBI GI.
/// It returns all matches found, each represented as a [`SequenceID`] variant, stored in
/// [`SeqIdList`]
///
/// The identifiers are sorted by database priority:
/// PDB > SwissProt > UniProtID > UniRef > RefSeq > GenBank > Ensembl > NCBI GI > NCBI Taxonomy.
///
/// # Examples
///
/// ```
/// use bioshell_seq::sequence::{parse_sequence_id, SeqId};
///
/// let ids = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9 2' cyclic ADP-D-ribose synthase [taxid=1310613]");
/// assert_eq!(ids.len(), 3);
/// assert!(matches!(ids[0], SeqId::UniProtID(_)));
/// assert!(matches!(ids[1], SeqId::TrEmbl(_)));
/// assert!(matches!(ids[2], SeqId::TaxId(_)));
/// assert_eq!(&ids.to_string(), "UniProtID|ABTIR_ACIB9|tr|A0A009IHW8|[taxid=1310613]");
/// let ids = parse_sequence_id(">ref|XP_001234567.1| hypothetical protein [Homo sapiens]");
/// assert_eq!(ids.len(), 1);
/// assert!(matches!(ids[0], SeqId::RefSeq(_)));
///
/// let ids = parse_sequence_id("NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome");
/// assert_eq!(&ids.to_string(), "RefSeq|NC_000913.3");
/// ```
pub fn parse_sequence_id(description: &str) -> SeqIdList {
    let patterns: &[(&str, fn(String) -> SeqId)] = &[
        (r"(?:pdb|\s+|\|)([0-9][A-Za-z0-9]{3})(?::[_]?[A-Za-z0-9]{0,3})?[ |]", |s| SeqId::PDB(s)),
        (r"\b([A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-\d+)?)\b", |s| SeqId::SwissProt(s)),
        (r"\b(UniRef\d{2,3}_[A-Z0-9]+)\b", |s| SeqId::UniRef(s)),
        (r"\b((?:NP|XP|WP|YP|XM|XR|NM|NR|NC)_[0-9]+\.\d+)\b", |s| SeqId::RefSeq(s)),
        (r"\bGI:(\d+)\b", |s| SeqId::NCBIGI(s)),
        (r"\bgi\|(\d+)\b", |s| SeqId::NCBIGI(s)),
        (r"\b(ENS[TPGR][0-9]{11})\b", |s| SeqId::Ensembl(s)),
        (r"\b([A-Z0-9]{3,6}_[A-Z0-9]{2,6})\b", |s| SeqId::UniProtID(s)),
        (r"\b([A-NR-Z0-9]{10})\b", |s| SeqId::TrEmbl(s)),
        (r"(?i:\[?taxid=(\d+))", |s| SeqId::TaxId(s)),
        (r"OX=(\d+)", |s| SeqId::TaxId(s)),
        (r"(?:\b|\|)gb\|([A-Z]{1,3}[0-9]{4,8}(?:\.\d+)?)\b", |s| SeqId::GenBank(s)),
    ];

    let mut found = Vec::new();
    let mut buffer = description.as_bytes().to_vec();

    for (pattern, constructor) in patterns {
        let re = Regex::new(pattern).unwrap();
        let desc_view = std::str::from_utf8(&buffer).unwrap();
        if let Some(caps) = re.captures(desc_view) {
            if let Some(m) = caps.get(1) {
                let matched_str = m.as_str().to_string();
                found.push(constructor(matched_str));
                // Mask out the matched portion with whitespace of same length
                for i in m.start()..m.end() {
                    buffer[i] = b' ';
                }
            }
        }
    }

    if found.is_empty() {
        found.push(SeqId::Default(description.split_whitespace().next().unwrap_or("").to_string()));
    }

    found.sort(); // Sorts by priority
    // eprintln!("Found: {:?}", &found);
    SeqIdList::from(found)
}

/// A typed wrapper around a list of sequence identifiers.
///
/// This struct provides ergonomic methods and traits for handling collections of `SeqId`,
/// including sorting by biological database priority and formatting into standard headers.
///
/// The formatting (`Display` / `.to_string()`) produces a `|`-separated header string such as:
/// `"PDB|1HHP:A|SwissProt|Q9NQX5|RefSeq|XP_123456.1"`
///
/// # Example
///
/// ```
/// use bioshell_seq::sequence::{SeqId, SeqIdList};
///
/// let mut ids = SeqIdList::from(vec![
///     SeqId::RefSeq("XP_123456.1".to_string()),
///     SeqId::SwissProt("Q9NQX5".to_string()),
///     SeqId::PDB("1HHP:A".to_string()),
/// ]);
///
/// ids.sort();  // Sort by biological database importance
///
/// // Format as a standard header string
/// let mut header = ids.to_string();
/// assert_eq!(header, "pdb|1HHP:A|sp|Q9NQX5|RefSeq|XP_123456.1");
///
/// println!("Formatted Header: {}", ids); // uses Display
/// ```
///
/// You can also iterate over `SeqIdList` or access its inner `Vec<SeqId>` via deref.

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SeqIdList(pub Vec<SeqId>);

impl SeqIdList {
    /// Sorts the identifiers in-place by database priority (PDB > SwissProt > ...).
    pub fn sort(&mut self) { self.0.sort(); }

    /// Returns a sanitized, filesystem-safe string suitable for use as a file name.
    ///
    /// This uses the internal sequence IDs and joins them by underscores, replacing or removing
    /// problematic characters (`|`, `:`, `,`, ` `, `/`, etc.) to ensure compatibility with
    /// both Windows and Unix-like file systems.
    ///
    /// # Example
    ///
    /// ```
    /// use bioshell_seq::sequence::{SeqId, SeqIdList};
    /// let ids = SeqIdList::from(vec![
    ///     SeqId::RefSeq("XP_123456.1".to_string()),
    ///     SeqId::SwissProt("Q9NQX5".to_string()),
    /// ]);
    /// assert_eq!(ids.file_name(), "RefSeq_XP_123456.1_sp_Q9NQX5");
    /// ```
    pub fn file_name(&self) -> String {
        if self.0.is_empty() {
            return "sequence_ids".to_string();
        }

        let name = self.0.iter().flat_map(|id| {
            let (label, value) = match id {
                SeqId::PDB(s) => ("pdb", s),
                SeqId::SwissProt(s) => ("sp", s),
                SeqId::UniProtID(s) => ("UniProtID", s),
                SeqId::UniRef(s) => ("UniRef", s),
                SeqId::RefSeq(s) => ("RefSeq", s),
                SeqId::GenBank(s) => ("GenBank", s),
                SeqId::Ensembl(s) => ("Ensembl", s),
                SeqId::NCBIGI(s) => ("NCBIGI", s),
                SeqId::TrEmbl(s) => ("tr", s),
                SeqId::TaxId(s) => ("taxid", s),
                SeqId::Default(s) => ("", s),
            };
            vec![label.to_string(), value.to_string()]
        })
            .collect::<Vec<_>>()
            .join("_");
        let name = name.trim_matches(&['_']);
        sanitize_filename(&name)
    }
}

impl From<Vec<SeqId>> for SeqIdList {
    fn from(vec: Vec<SeqId>) -> Self { SeqIdList(vec) }
}

impl IntoIterator for SeqIdList {
    type Item = SeqId;
    type IntoIter = std::vec::IntoIter<SeqId>;

    fn into_iter(self) -> Self::IntoIter { self.0.into_iter() }
}

impl<'a> IntoIterator for &'a SeqIdList {
    type Item = &'a SeqId;
    type IntoIter = std::slice::Iter<'a, SeqId>;

    fn into_iter(self) -> Self::IntoIter { self.0.iter() }
}

impl Deref for SeqIdList {
    type Target = Vec<SeqId>;

    fn deref(&self) -> &Self::Target { &self.0 }
}

impl DerefMut for SeqIdList {
    fn deref_mut(&mut self) -> &mut Self::Target { &mut self.0 }
}


impl fmt::Display for SeqIdList {
    /// Allows formatting the entire list using `.to_string()` or `{}`.
    ///
    /// The output is a `|`-separated string of alternating database names and identifiers,
    /// e.g., `PDB|1HHP:A|sp|P12345|RefSeq|XP_001234567.1`.
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::sequence::{SeqId, SeqIdList};
    ///
    /// let ids = vec![
    ///     SeqId::SwissProt("Q9NQX5".to_string()),
    ///     SeqId::RefSeq("XP_001234567.1".to_string()),
    ///     SeqId::PDB("1HHP:A".to_string()),
    /// ];
    ///
    /// let mut ids = SeqIdList::from(ids);
    /// ids.sort();
    /// let header = ids.to_string();
    /// assert_eq!(header, "pdb|1HHP:A|sp|Q9NQX5|RefSeq|XP_001234567.1");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, id) in self.0.iter().enumerate() {
            if i > 0 { write!(f, "|")?; }
            write!(f, "{id}")?;
        }
        Ok(())
    }
}
