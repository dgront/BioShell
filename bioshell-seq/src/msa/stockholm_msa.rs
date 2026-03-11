
use std::collections::HashMap;
use std::fmt;
use std::ops::Deref;
use crate::msa::MSA;
use crate::sequence::Sequence;

/// [`StockholmMSA`] is an [`MSA`] struct with additional annotations.
///
/// The struct stores a regular [`MSA`] together with Stockholm annotations parsed from a file in the Stockholm format, e.g.:
///
/// ```text
#[doc = include_str!("../../tests/test_files/4Fe-4S-example.sto")]
/// ```
///
/// The first line of a Stockholm file (``.sto``, or ``.stk``) states the format and version identifier,
/// currently ``# STOCKHOLM 1.0``. The header is followed by mark-up lines beginning with ``#``.
/// These mark-up lines can annotate features of the alignment file (``#=GF``, generic per-file annotation),
/// or features of the aligned sequences (``#=GS``, generic per-sequence annotation).
/// The sequence alignment itself is a series of lines with sequence names (typically in the form name/start-end)
/// followed by a space and the aligned sequence. A line with two forward slashes (``//``) indicates the end of the alignment.
///
/// For the detailed description and a list of allowed annotations, see [Stockholm specification](https://sonnhammer.sbc.su.se/Stockholm.html)
///
/// Having the file loaded, one can access both the aligned sequences and the respective annotations:
///
/// ```rust
/// use bioshell_seq::msa::StockholmMSA;
/// # use bioshell_seq::SequenceError;
/// # fn main() -> Result<(), SequenceError> {
/// # use std::fs::File;
/// # use std::io::BufReader;
///
/// let file = File::open("tests/test_files/4Fe-4S-example.sto")?;
/// let mut reader = BufReader::new(file);
///
/// let sto = StockholmMSA::from_stockholm_reader(&mut reader)?;
/// assert_eq!(sto.n_seq(), 6);
/// assert!(sto.len() > 0);
///
/// # Ok(())
/// # }
/// ```
/// The wrapped [`MSA`] is exposed by delegation via `Deref`, so this type behaves like [`MSA`].
///
#[derive(Debug, Clone)]
pub struct StockholmMSA {
    pub(crate) msa: MSA,
    pub(crate) gf: HashMap<GfEntryType, Vec<String>>,
    pub(crate) gs: HashMap<String, HashMap<GsEntryType, Vec<String>>>,
    pub(crate) gc: HashMap<GcEntryType, String>,
    pub(crate) gr: HashMap<String, HashMap<GrEntryType, String>>,
}

impl StockholmMSA {

    /// Returns the wrapped `MSA`.
    pub fn msa(&self) -> &MSA { &self.msa }

    /// Consumes `self` and returns the wrapped `MSA`.
    pub fn into_msa(self) -> MSA {
        self.msa
    }

    /// Returns all values stored for a given GF key.
    pub fn gf(&self, key: GfEntryType) -> Option<&[String]> {
        self.gf.get(&key).map(Vec::as_slice)
    }

    /// Returns all values stored for a given GS key of a given sequence.
    pub fn gs(&self, seq_id: &str, key: GsEntryType) -> Option<&[String]> {
        self.gs.get(seq_id)?.get(&key).map(Vec::as_slice)
    }

    /// Returns the GC value for a given key.
    pub fn gc(&self, key: GcEntryType) -> Option<&str> {
        self.gc.get(&key).map(String::as_str)
    }

    /// Returns the GR value for a given sequence and key.
    pub fn gr(&self, seq_id: &str, key: GrEntryType) -> Option<&str> {
        self.gr.get(seq_id)?.get(&key).map(String::as_str)
    }

    /// Alignment length, i.e. the length of each of the aligned sequence including gaps
    pub fn len(&self) -> usize {
        self.msa.len()
    }

    /// Returns true if the alignment is empty.
    pub fn is_empty(&self) -> bool {
        self.n_seq() == 0 || self.len() == 0
    }

    /// Number of sequences aligned.
    pub fn n_seq(&self) -> usize {
        self.msa.n_seq()
    }

    /// Returns immutable access to aligned sequences.
    pub fn sequences(&self) -> &Vec<Sequence> {
        self.msa.sequences()
    }
}

impl Deref for StockholmMSA {
    type Target = MSA;

    fn deref(&self) -> &Self::Target { &self.msa }
}

impl From<MSA> for StockholmMSA {
    /// Allows to create a [`StockholmMSA`] struct from [`MSA`]
    ///
    fn from(msa: MSA) -> Self {
        Self {
            msa,
            gf: HashMap::new(),
            gs: HashMap::new(),
            gc: HashMap::new(),
            gr: HashMap::new(),
        }
    }
}

/// Generate a Stockholm entry-type enum together with parsing/formatting helpers.
macro_rules! stockholm_entry_type {
    (
        $(#[$enum_meta:meta])*
        $vis:vis enum $name:ident {
            $(
                $(#[$variant_meta:meta])*
                $variant:ident => $text:literal
            ),+ $(,)?
        }
    ) => {
        $(#[$enum_meta])*
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        $vis enum $name {
            $(
                $(#[$variant_meta])*
                $variant,
            )+
            /// Any non-standard annotation key.
            Other(String),
        }

        impl $name {
            /// Returns the canonical Stockholm label.
            pub fn as_str(&self) -> &str {
                match self {
                    $(Self::$variant => $text,)+
                    Self::Other(s) => s.as_str(),
                }
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(self.as_str())
            }
        }

        impl FromStr for $name {
            type Err = core::convert::Infallible;

            fn from_str(s: &str) -> Result<Self, Self::Err> {
                Ok(match s {
                    $($text => Self::$variant,)+
                    other => Self::Other(other.to_string()),
                })
            }
        }

        impl From<&str> for $name {
            fn from(value: &str) -> Self {
                value.parse().unwrap()
            }
        }
    };
}

use std::str::FromStr;

stockholm_entry_type! {
    /// Stockholm `GF` (global/file) annotation entry types.
    pub enum GfEntryType {
        /// Alignment identifier.
        ID => "ID",
        /// Accession.
        AC => "AC",
        /// Description.
        DE => "DE",
        /// Author line.
        AU => "AU",
        /// Seed/source information.
        SE => "SE",
        /// Alignment source/method note.
        SS => "SS",
        /// Gathering cutoff.
        GA => "GA",
        /// Trusted cutoff.
        TC => "TC",
        /// Noise cutoff.
        NC => "NC",
        /// Build method.
        BM => "BM",
        /// Search method.
        SM => "SM",
        /// Entry type.
        TP => "TP",
        /// Database reference.
        DR => "DR",
        /// Comment.
        CC => "CC",
        /// Reference number.
        RN => "RN",
        /// Reference Medline/PubMed id.
        RM => "RM",
        /// Reference title.
        RT => "RT",
        /// Reference authors.
        RA => "RA",
        /// Reference location.
        RL => "RL",
        /// Reference comment.
        RC => "RC",
        /// Keywords.
        KW => "KW",
    }
}

stockholm_entry_type! {
    /// Stockholm `GS` (per-sequence) annotation entry types.
    pub enum GsEntryType {
        /// Sequence accession.
        AC => "AC",
        /// Sequence description.
        DE => "DE",
        /// Database reference.
        DR => "DR",
        /// Organism name.
        OS => "OS",
        /// Organism classification.
        OC => "OC",
        /// Sequence location.
        LO => "LO",
    }
}

stockholm_entry_type! {
    /// Stockholm `GC` (per-column) annotation entry types.
    pub enum GcEntryType {
        /// Consensus secondary structure.
        SSCons => "SS_cons",
        /// Reference annotation.
        RF => "RF",
        /// Model mask.
        MM => "MM",
    }
}

stockholm_entry_type! {
    /// Stockholm `GR` (per-residue, per-sequence) annotation entry types.
    pub enum GrEntryType {
        /// Residue secondary structure.
        SS => "SS",
        /// Surface accessibility.
        SA => "SA",
        /// Transmembrane annotation.
        TM => "TM",
        /// Posterior probability.
        PP => "PP",
        /// Ligand interaction.
        LI => "LI",
        /// Active site.
        AS => "AS",
        /// Intron annotation.
        IN => "IN",
    }
}

