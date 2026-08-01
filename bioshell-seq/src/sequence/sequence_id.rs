use std::cmp::Ordering;
use std::fmt;
use std::ops::{Deref, DerefMut};
use std::sync::LazyLock;
use regex::{Captures, Regex};

use bioshell_core::io::sanitize_filename;

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
/// The [`parse_sequence_id()`](parse_sequence_id) function can be used to parse a sequence description string into a list of `SeqId` variants.
/// For convenience,  [`parse_sequence_id()`](parse_sequence_id) returns [`SeqIdList`].
/// ```
/// use bioshell_seq::sequence::{parse_sequence_id, SeqId};
/// let ids = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9 2' cyclic ADP-D-ribose synthase [taxid=1310613]");
/// assert_eq!(ids.len(), 3);
/// assert!(matches!(ids[0], SeqId::SwissProt(_)));
/// assert!(matches!(ids[1], SeqId::UniProtEntry(_)));
/// assert!(matches!(ids[2], SeqId::TaxId(_)));
/// assert_eq!(&ids.to_string(), "sp|A0A009IHW8|ABTIR_ACIB9|taxid=1310613");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SeqId {
    /// CYP-id - the identifier for a P450 cytochrome used by the P450Atlas
    CypId(String),

    /// PDB (Protein Data Bank) entry ID, optionally with chain (e.g., "1HHP", "1HHP:A", "pdb_00002gb1").
    ///
    /// This variant captures both the new 12-characters long identifiers
    /// and the legacy type, which are 4-character alphanumeric codes. An optional chain identifier can be included after a colon.
    PDB(String),

    /// Accession ID to SwissProt part of UniProtKB (e.g., "sp|P12345", "sp|Q9NQX5-2").
    SwissProt(String),

    /// TrEmbl accession ID  of UniProtKB (e.g., "tr|A0A009IHW8|" or "tr|P12345").
    TrEmbl(String),

    /// UniProtKB accession ID when it's not specified either it's SwissProt or TrEMBL.
    ///
    /// An accession number (AC) is assigned to each sequence upon inclusion into UniProtKB.
    /// Accession numbers are stable from release to release. For details, see the official
    /// specification:
    /// [UniProt accession number format](https://www.uniprot.org/help/accession_numbers).
    UniProtKB(String),

    /// UniParc accession, e.g. UPI001E6762CD
    ///
    /// The UniParc identifier (UPI) is the unique identifier assigned to a distinct
    /// protein sequence. It consists of the characters “UPI” followed
    /// by 10 hexadecimal characters (0–9, A–F). A UPI is stable across releases
    /// [The UniParc identifier (UPI)](https://www.uniprot.org/help/uniparc_id).
    UniParc(String),

    /// UniProtKB entry name (e.g., "002L_FRG3G").
    ///
    /// The 'Entry name' is a unique identifier, often containing biologically relevant information.
    /// For details, see the official page:
    /// [difference between an accession number (AC) and the entry name](https://www.uniprot.org/help/difference_accession_entryname).
    UniProtEntry(String),

    /// UniRef (UniProt Reference Clusters), e.g., "UniRef100_P12345".
    UniRef(String),

    /// RefSeq protein or transcript accession from NCBI (e.g., "XP_123456.1", "NM_001256789.2").
    ///
    /// The RefSeq ID string starts with a two-letter prefix indicating the type of sequence, e.g., "XP" for predicted protein
    /// The list of valid prefixes can be found in the
    /// [NCBI RefSeq accession prefix table](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/).
    RefSeq(String),

    /// GenBank or EMBL accession (e.g., "AB123456.1", "U49845").
    ///
    /// Specification can be found on the NCBI website:
    /// [Accession Number prefixes](https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/)
    GenBank(String),

    /// DNA Data Bank of Japan entry, e.g. `dbj|BAE93469.1`
    DDBJ(String),

    /// Ensembl gene, transcript, or protein ID (e.g., "ENSG00000139618", "ENST00000331789").
    Ensembl(String),

    /// NCBI GI number, now obsolete ("gi|12345678" or "GI:12345678").
    NCBIGI(String),

    /// KEGG database id.
    KEGG(String),

    /// NCBI Taxonomy ID (e.g., "taxid=9606", "[taxid=9606]", "OX=9606").
    ///
    /// This variant captures the taxonomy annotation both in the NCBI and the UniProt format,
    /// i.e. "taxid=9606" and "OX=9606", respectively.
    TaxId(String),

    /// Organism scientific name (e.g., "Homo sapiens")
    ///
    /// This variant captures both the species annotation in the UniProt and NCBI formats,
    /// e.g. "OS=9606" or "[Homo sapiens]"
    Organism(String),

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
            SeqId::CypId(_) => 0,
            SeqId::PDB(_) => 1,
            SeqId::SwissProt(_) => 2,
            SeqId::UniProtKB(_) => 3,
            SeqId::TrEmbl(_) => 4,
            SeqId::UniProtEntry(_) => 5,
            SeqId::UniParc(_) => 6,
            SeqId::UniRef(_) => 7,
            SeqId::RefSeq(_) => 8,
            SeqId::GenBank(_) => 9,
            SeqId::Ensembl(_) => 10,
            SeqId::DDBJ(_) => 11,
            SeqId::NCBIGI(_) => 12,
            SeqId::KEGG(_) => 13,
            SeqId::Default(_) => 14,
            SeqId::TaxId(_) => 15,
            SeqId::Organism(_) => 16,
        }
    }
    
    pub fn value(&self) -> &str {
        match self {
            SeqId::CypId(v) => v,
            SeqId::PDB(v) => v,
            SeqId::SwissProt(v) => v,
            SeqId::UniProtKB(v) => v,
            SeqId::TrEmbl(v) => v,
            SeqId::UniProtEntry(v) => v,
            SeqId::UniParc(v) => v,
            SeqId::UniRef(v) => v,
            SeqId::RefSeq(v) => v,
            SeqId::GenBank(v) => v,
            SeqId::Ensembl(v) => v,
            SeqId::DDBJ(v) => v,
            SeqId::NCBIGI(v) => v,
            SeqId::Default(v) => v,
            SeqId::TaxId(v) => v,
            SeqId::Organism(v) => v,
            SeqId::KEGG(v) => v,
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
    /// assert_eq!(header, "ref|XP_001234567.1");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SeqId::PDB(s) if s.starts_with("pdb_") => write!(f, "{}", s),
            SeqId::PDB(s) => write!(f, "pdb|{}", s),
            SeqId::SwissProt(s) => write!(f, "sp|{}", s),
            SeqId::TrEmbl(s) => write!(f, "tr|{}", s),
            SeqId::UniProtKB(s) => write!(f, "UniProt|{}", s),
            SeqId::UniParc(s) => write!(f, "UniParc|{}", s),
            SeqId::UniProtEntry(s) => write!(f, "{}", s),
            SeqId::UniRef(s) => write!(f, "{}", s),
            SeqId::KEGG(s) => write!(f, "{}", s),
            SeqId::RefSeq(s) => write!(f, "ref|{}", s),
            SeqId::GenBank(s) => write!(f, "gb|{}", s),
            SeqId::Ensembl(s) => write!(f, "Ensembl|{}", s),
            SeqId::NCBIGI(s) => write!(f, "gi|{}", s),
            SeqId::TaxId(s) => write!(f, "taxid={}", s),
            SeqId::Default(s) => write!(f, "{}", s),
            SeqId::Organism(s) => write!(f, "[organism={}]", s),
            SeqId::CypId(s) => write!(f, "{}", s),
            SeqId::DDBJ(s) => write!(f, "dbj|{}", s),
        }
    }
}

/// Extract sequence identifiers from a free-text description string.
///
/// This function scans the input for standard database identifiers such as those from
/// PDB, SwissProt/TrEMBL, UniRef, RefSeq, GenBank/EMBL, Ensembl, and NCBI GI.
/// It returns all matches found, each represented as a [`SeqId`] variant, stored in
/// [`SeqIdList`]
///
/// The returned identifiers are stored in the order as they appeared in the given description string.
/// One can use [`sort()`](SeqIdList::sort) method to sort the entries by the database priority:
/// PDB > SwissProt > UniProtID > UniRef > RefSeq > GenBank > Ensembl > NCBI GI > NCBI Taxonomy.
///
/// # Examples
///
/// ```
/// use bioshell_seq::sequence::{parse_sequence_id, SeqId};
///
/// let ids = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9 2' cyclic ADP-D-ribose synthase [taxid=1310613]");
/// assert_eq!(ids.len(), 3);
/// assert!(matches!(ids[0], SeqId::SwissProt(_)));
/// assert!(matches!(ids[1], SeqId::UniProtEntry(_)));
/// assert!(matches!(ids[2], SeqId::TaxId(_)));
/// assert_eq!(&ids.to_string(), "sp|A0A009IHW8|ABTIR_ACIB9|taxid=1310613");
/// let ids = parse_sequence_id(">ref|XP_001234567.1| hypothetical protein [Homo sapiens]");
/// assert_eq!(ids.len(), 2);
/// assert!(matches!(ids[0], SeqId::RefSeq(_)));
/// assert!(matches!(ids[1], SeqId::Organism(_)));
///
/// let ids = parse_sequence_id("NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome");
/// assert_eq!(&ids.to_string(), "ref|NC_000913.3");
/// ```
pub fn parse_sequence_id(description: &str) -> SeqIdList {

    let desc = expand_taxids(description);
    let mut found_with_positions = Vec::new();
    let mut buffer = desc.as_bytes().to_vec();

    for (pattern, constructor) in SEQUENCE_ID_PATTERNS.iter() {
        let matches: Vec<(usize, usize, String)> = {
            let desc_view = std::str::from_utf8(&buffer).unwrap();

            pattern
                .captures_iter(desc_view)
                .filter_map(|caps| {
                    let whole = caps.get(0)?;
                    let id = caps.get(1)?;

                    Some((
                        whole.start(),
                        whole.end(),
                        id.as_str().to_owned(),
                    ))
                })
                .collect()
        };

        for (start, end, matched_str) in matches {
            found_with_positions.push((start, constructor(matched_str)));
            buffer[start..end].fill(b' ');
        }
    }

    found_with_positions.sort_by_key(|(start, _)| *start);

    let mut found = found_with_positions
        .into_iter()
        .map(|(_, id)| id)
        .collect::<Vec<_>>();

    if found.is_empty() {
        found.push(SeqId::Default(description.split_whitespace().next().unwrap_or("").to_string()));
    }

    // found.sort(); // Sorts by priority
    SeqIdList::from(found)
}

static TAXID_LIST_RE: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"(?i)\btaxid(?:=|\|)(\d+(?:,\d+)*)(?:\|)?[ \t]*")
        .expect("invalid TAXID_LIST_RE")
});

fn expand_taxids(text: &str) -> String {

    if !text.contains("taxid") {
        return text.to_owned();
    }

    TAXID_LIST_RE.replace_all(text, |caps: &Captures| {
        caps[1]
            .split(',')
            .map(|id| format!("taxid|{id}|"))
            .collect::<Vec<_>>()
            .join(" ") + " "
    }).trim_end().to_string()
}


/// A list of sequence identifiers.
///
/// This struct provides ergonomic methods to operate on `SeqId`, such as:
///  - sorting by database priority,
///  - formatting into a standard header string,
///  - generating a filesystem-safe file name.
///
/// The formatting (`Display` / `.to_string()`) produces a `|`-separated header string such as:
/// `"PDB|1HHP:A|sp|Q9NQX5|ref|XP_123456.1"`
///
/// # Examples
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
/// assert_eq!(header, "pdb|1HHP:A|sp|Q9NQX5|ref|XP_123456.1");
///
/// println!("Formatted Header: {}", ids); // uses Display
/// ```
///
/// A [`SeqIdList`] object is typically obtained by calling [`parse_sequence_id()`]:
/// ```
/// # use bioshell_seq::sequence::parse_sequence_id;
/// let id_list = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9|taxid=1310613");
/// ```
///
/// You can also iterate over [`SeqIdList`] or access its inner `Vec<SeqId>` via deref.
/// E.g. to check if it contains particular database id (PDB in this example), say:
/// ```
/// # use bioshell_seq::sequence::{parse_sequence_id, SeqId};
/// let ids = parse_sequence_id("sp|A0A009IHW8|ABTIR_ACIB9|taxid=1310613");
/// let perhaps_pdb = ids.iter().find_map(|id| match id {
///         SeqId::PDB(value) => Some(value.as_str()),
///         _ => None,
///     });
/// # assert!(perhaps_pdb.is_none());
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SeqIdList(pub Vec<SeqId>);

impl SeqIdList {
    /// Sorts the identifiers in-place by decreasing database priority.
    ///
    /// The pre-defined order is: PDB > SwissProt > UniProtID > UniRef > RefSeq > GenBank > Ensembl > TrEmbl > NCBI GI > NCBI Taxonomy;
    /// i.e., the most important identifiers (e.g. PDB id) appear first in the list.
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
    /// assert_eq!(ids.file_name(), "ref_XP_123456.1_sp_Q9NQX5");
    /// ```
    pub fn file_name(&self) -> String {
        if self.0.is_empty() {
            return "sequence_ids".to_string();
        }

        let name = format!("{}", self)
            .replace('|', "_")
            .replace(']', "")
            .replace(' ', "_")
            .replace("[organism=", "");
        let name = name.trim_matches(&['_','|']);
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
    /// assert_eq!(header, "pdb|1HHP:A|sp|Q9NQX5|ref|XP_001234567.1");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, id) in self.0.iter().enumerate() {
            if i > 0 {
                write!(f, "{}", if matches!(id, SeqId::Organism(_)) { " " } else { "|" })?;
            }
            write!(f, "{id}")?;
        }
        Ok(())
    }
}

type SeqIdConstructor = fn(String) -> SeqId;
type CompiledPattern = (Regex, SeqIdConstructor);

macro_rules! regex_patterns {
    ($($pattern:expr => $constructor:expr),* $(,)?) => {
        LazyLock::new(|| {
            vec![
                $(
                    (
                        Regex::new($pattern).unwrap_or_else(|error| {
                            panic!(
                                "invalid sequence ID regex {:?}: {}",
                                $pattern,
                                error
                            )
                        }),
                        $constructor as SeqIdConstructor,
                    )
                ),*
            ]
        })
    };
}

static SEQUENCE_ID_PATTERNS: LazyLock<Vec<CompiledPattern>> = regex_patterns![

            // tax-id can be in various formats, e.g., [taxid=9606], taxid=9606, taxid|9606, TaxID=9606, OX=9606
            r"(?i:\[?taxid=(\d+))" => |s| SeqId::TaxId(s),
            r"(?i:\[?TaxID=(\d+))" => |s| SeqId::TaxId(s),
            r"OX=(\d+)" => |s| SeqId::TaxId(s),
            r"(?i:(?:\b|\|)taxid\|(\d+))" => |s| SeqId::TaxId(s),
            // KEGG id
            r"(?:^|[|>]|\b)([a-z][a-z0-9]{2,4}:[A-Za-z0-9_.-]+)(?:$|[|>]|\b)" => |s| SeqId::KEGG(s),

            r"(?:^|pdb|\s+|\|)([0-9][A-Za-z0-9]{3}(?::[_]?[A-Za-z0-9]{1,3})?)(?:$|[ |])" => |s| SeqId::PDB(s),
            r"\b(pdb_[A-Za-z0-9]{8}(?::[_]?[A-Za-z0-9]{1,3})?)(?:$|[ |])" => |s| SeqId::PDB(s),

            // RefSeq must be checked early because it's similar to UniProt entry name
            r"(?:\b|\|\>|_)((?:AC|NC|NG|NT|NW|NZ|NM|NR|XM|XR|AP|NP|YP|XP|WP)_[0-9]+\.\d+)\b" => |s| SeqId::RefSeq(s),

            // SwissProt accession with explicit prefix, all variants
            r"(?:\b|\|\>)sp[|.]([A-Z0-9]{6}|[A-Z0-9]{10})(?:-\d+)?[|.]" => |s| SeqId::SwissProt(s),
            // trEmbl accession with explicit prefix, all variants
            r"(?:\b|\|)tr[|.]([A-Z0-9]{6}|[A-Z0-9]{10})(?:-\d+)?[|.]" => |s| SeqId::TrEmbl(s),
            // UniProt entry name, e.g., "002L_FRG3G"
            r"(?:\b|\|)([A-Z0-9]{3,}_[A-Z0-9]{3,5})\b" => |s| SeqId::UniProtEntry(s),
            // UniProtKB accession: legacy 6-character versions and the new 10-character form
            r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\b" => |s| SeqId::UniProtKB(s),
            // UniRef
            r"\b(UniRef\d{2,3}_[A-Z0-9]+)\b" => |s| SeqId::UniRef(s),
            // UniParc ID consists of the characters “UPI” followed by 10 hexadecimal characters (0–9, A–F).
            r"\b(UPI[0-9A-F]{10})\b" => |s| SeqId::UniParc(s),

            r"\bGI:(\d+)\b" => |s| SeqId::NCBIGI(s),
            r"\bgi\|(\d+)\b" => |s| SeqId::NCBIGI(s),
            r"\b(ENS[TPGR][0-9]{11})\b" => |s| SeqId::Ensembl(s),

            // INSDC / DDBJ explicitly prefixed with "dbj|", e.g., "dbj|BAE93469.1"
            r"(?:\b|\|)dbj\|([A-Z]{3}[0-9]{5}(?:\.\d+)?)\b" => |s| SeqId::DDBJ(s),
            // INSDC / GenBank explicitly prefixed with "gb|", e.g., "gb|AB123456.1"
            r"(?:\b|\||_)gb\|([A-Z]{1,3}[0-9]{4,8}(?:\.\d+)?)\b" => |s| SeqId::GenBank(s),
            // INSDC / GenBank nucleotide: 1 letter + 5 digits
            r"(?:\b|\||_)([A-Z][0-9]{5}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),
            // INSDC / GenBank nucleotide: 2 letters + 6 digits
            r"(?:\b|\||_)([A-Z]{2}[0-9]{6}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),
            // INSDC / GenBank nucleotide: 2 letters + 8 digits
            r"(?:\b|\||_)([A-Z]{2}[0-9]{8}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),
            // INSDC / GenBank protein: 3 letters + 5 digits
            r"(?:\b|\||_)([A-Z]{3}[0-9]{5}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),
            // INSDC / GenBank protein: 3 letters + 7 digits
            r"(?:\b|\||_)([A-Z]{3}[0-9]{7}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),
            // INSDC/GenBank WGS-style nucleotide accession
            r"(?:\b|\||_)([A-Z]{4}[0-9]{8,10}(?:\.\d+)?)(?:$|[^\w]|_)" => |s| SeqId::GenBank(s),

            // organism scientific name in square brackets
            r"\[organism=([[:alpha:]][[:alnum:]. ]*)\]" => |s| SeqId::Organism(s),
            r"\[([[:alpha:]][[:alnum:]. ]*)\]" => |s| SeqId::Organism(s),
            // CYP-id - two variants: lowercase and uppercase
            r"(?:^|[^\w]|_)(CYP[0-9]+[A-Z]{1,3}[0-9]+[a-z]?(?:v[0-9]{1,2})?(?:P(?:[0-9]+|[NC])?)?X?(?:_[A-Z]{1,4})?)(?:$|[^\w]|_)" => |s| SeqId::CypId(s),
            r"(?:^|[^\w]|_)(Cyp[0-9]+[a-z]{1,3}[0-9]+[a-z]?(?:v[0-9]{1,2})?(?:P(?:[0-9]+|[NC])?)?X?(?:_[A-Z]{1,4})?)(?:$|[^\w]|_)" => |s| SeqId::CypId(s),
        ];
