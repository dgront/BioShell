use regex::{Error, Regex};

use crate::sequence::{count_residue_type, Sequence};
use crate::SequenceError;

/// Says `true` or `false` for a [`Sequence`](Sequence) object
///
/// Sequence filters are be used to select sequences from a larger pool
pub trait SequenceFilter {
    /// Returns `true` or `false` for a given `sequence` object
    fn filter(&self, sequence: &Sequence) -> bool;
}

/// Returns `true` if a description of a given [`Sequence`](Sequence) contains a given substring
pub struct DescriptionContains { pub substring: String}

impl SequenceFilter for DescriptionContains {
    fn filter(&self, sequence: &Sequence) -> bool {
        sequence.description().contains(self.substring.as_str())
    }
}

/// Always returns `true`; serves as the defafult filter that pass every sequence
pub struct AlwaysTrue;

impl SequenceFilter for AlwaysTrue {
    fn filter(&self, _: &Sequence) -> bool { true }
}


/// Reverse the result of a given [`SequenceFilter`](SequenceFilter)
///
/// ```
/// use bioshell_seq::sequence::{AlwaysTrue, LogicalNot, Sequence, SequenceFilter};
/// let f = LogicalNot::new(AlwaysTrue);
/// let sequence = Sequence::from_str("test_seq", "MRAA");
/// assert!(!f.filter(&sequence));
/// ```
pub struct LogicalNot<F> { pub f: F}

impl<F: SequenceFilter> SequenceFilter for LogicalNot<F> {
    fn filter(&self, sequence: &Sequence) -> bool { !self.f.filter(sequence) }
}

/// Returns `true` if the length of a given [`Sequence`](Sequence) is within a certain range
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, SequenceLengthWithinRange, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAGIA");
/// let filter = SequenceLengthWithinRange{from: 3,to: 5};
/// assert_eq!(filter.filter(&sequence1), true);
/// assert_eq!(filter.filter(&sequence2), false);
/// ```
pub struct SequenceLengthWithinRange {
    pub from: usize,
    pub to: usize
}

impl SequenceFilter for SequenceLengthWithinRange {
    fn filter(&self, sequence: &Sequence) -> bool {
        sequence.len() >= self.from && sequence.len() <= self.to
    }
}

/// Returns `true` if a given [`Sequence`](Sequence) contains at least `n` unknown residues, marked as `'X'`
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, ContainsX, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAXXA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAXXXA");
/// let filter = ContainsX{min_x: 3};
/// assert_eq!(filter.filter(&sequence1),false);
/// assert_eq!(filter.filter(&sequence2),true);
/// ```
pub struct ContainsX {
    pub min_x: usize,
}

impl SequenceFilter for ContainsX {
    fn filter(&self, sequence: &Sequence) -> bool {
        count_residue_type(sequence, 'X') >= self.min_x
    }
}

/// Returns `true` if a given [`Sequence`](Sequence) contains at least `f` fraction of unknown residues, marked as `'X'`
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, FractionX, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAXXA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAGSXA");
/// let filter = FractionX{min_x_fraction: 0.25};
/// assert_eq!(filter.filter(&sequence1), true);
/// assert_eq!(filter.filter(&sequence2), false);
/// ```
pub struct FractionX {
    pub min_x_fraction: f64,
}

impl SequenceFilter for FractionX {
    fn filter(&self, sequence: &Sequence) -> bool {
        count_residue_type(sequence, 'X') as f64 >= self.min_x_fraction * sequence.len() as f64
    }
}

/// Returns `true` if a given [`Sequence`](Sequence) is a valid DNA or RNA sequence
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, IsNucleic, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "CGCGTATACGCG");
/// let sequence2 = Sequence::from_str("test_seq", "cgcgtatacgc");
/// let sequence3 = Sequence::from_str("test_seq", "CGATAGS");
/// let filter = IsNucleic;
/// assert!(filter.filter(&sequence1));
/// assert!(filter.filter(&sequence2));
/// assert!(!filter.filter(&sequence3));
/// ```
pub struct IsNucleic;

macro_rules! is_nucleotide {
    ($c: expr) => {
        match $c {
            b'A' | b'C' | b'T' | b'G' | b'U' | b'X' | b'a' | b'c' | b't' | b'g' | b'u' | b'x' => true,
            _ => false,
        }
    };
}

impl SequenceFilter for IsNucleic {
    fn filter(&self, sequence: &Sequence) -> bool {
        sequence.as_u8().iter().all(|l| is_nucleotide!(l))
    }
}

/// Returns `true` if a given [`Sequence`](Sequence) is a valid protein sequence
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, IsProtein, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAGSXA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAGSXA*");
/// let sequence3 = Sequence::from_str("test_seq", "MRAG!SXA");
/// let filter = IsProtein;
/// assert!(filter.filter(&sequence1));     // 'X' is a legal amino acid for IsProtein() filter
/// assert!(!filter.filter(&sequence2));    // ... but '*' is not allowed
/// assert!(!filter.filter(&sequence3));    // Neither does '!'
/// ```
pub struct IsProtein;
macro_rules! is_amino_acid {
    ($c: expr) => {
        match $c {
            b'A' | b'R' | b'N' | b'D' | b'C' | b'E' | b'Q' | b'G' | b'H' | b'I' | b'L' | b'K' | b'M' | b'F' | b'P' | b'S' | b'T' | b'W' | b'Y' | b'V' | b'X' => true,
            _ => false,
        }
    };
}

impl SequenceFilter for IsProtein {
    fn filter(&self, sequence: &Sequence) -> bool {
        sequence.as_u8().iter().all(|l| is_amino_acid!(l))
    }
}

/// Returns `true` if a given [`Sequence`](Sequence) contains no more than ``max_length`` residues
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, ShorterThan, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAGSXA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAGS");
/// let filter = ShorterThan{max_length: 6};
/// assert!(!filter.filter(&sequence1));
/// assert!(filter.filter(&sequence2));
/// ```
pub struct ShorterThan { pub max_length : usize }

impl SequenceFilter for ShorterThan {
    fn filter(&self, sequence: &Sequence) -> bool { sequence.len() <= self.max_length }
}

/// Returns `true` if a given [`Sequence`](Sequence) contains at least ``max_length`` residues
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, LongerThan, SequenceFilter};
/// let sequence1 = Sequence::from_str("test_seq", "MRAGSXA");
/// let sequence2 = Sequence::from_str("test_seq", "MRAGS");
/// let filter = LongerThan{min_length: 0};
/// assert!(filter.filter(&sequence1));
/// assert!(!filter.filter(&sequence2));
/// ```
pub struct LongerThan { pub min_length : usize }

impl SequenceFilter for LongerThan {
    fn filter(&self, sequence: &Sequence) -> bool { sequence.len() >= self.min_length }
}

/// Returns `true` if a given [`Sequence`](Sequence) contains a given sequence motif
///
/// # Examples
/// ```rust
/// use bioshell_seq::sequence::{Sequence, HasSequenceMotif, SequenceFilter};
/// use bioshell_seq::SequenceError;
/// # fn main() -> Result<(), SequenceError> {
/// let sequence1 = Sequence::from_str("no_motifs", "MRAGS");
/// let sequence2 = Sequence::from_str("p450_motifs", "PERFAGEILR");
/// let exxr_motif = HasSequenceMotif::new("ExxR")?;
/// assert!(!exxr_motif.filter(&sequence1));
/// assert!(exxr_motif.filter(&sequence2));
/// let perf_motif = HasSequenceMotif::new("PxR[FD]")?;
/// assert!(perf_motif.filter(&sequence2));
/// let perf_or_exxr = HasSequenceMotif::new("(ExxR|PxR[FD])")?;
/// assert!(perf_or_exxr.filter(&sequence2));
/// # Ok(())
/// # }
/// ```
pub struct HasSequenceMotif {  motif_regex : Regex }

impl HasSequenceMotif {
    pub fn new(pattern_str: &str) -> Result<HasSequenceMotif, SequenceError> {
        let regex_str = HasSequenceMotif::motif_to_regex(pattern_str);
        match Regex::new(&regex_str) {
            Ok(re_struct) => { Ok(HasSequenceMotif { motif_regex: re_struct }) }
            Err(_) => { Err(SequenceError::IncorrectSequencePattern { pattern: pattern_str.to_string() }) }
        }
    }

    /// Transforms a sequence motif pattern into a regex string.
    /// For example: "C-x(2,4)-C-[IL]-x(3)|(ExxR|PxR)" becomes "C.{2,4}C[IL].{3}|(E..R|P.R)".
    fn motif_to_regex(pattern: &str) -> String {
        // Regex pattern to find occurrences of x(n,m) or x(n)
        let re_repeats = Regex::new(r"x\((\d+),?(\d*)\)").unwrap();

        // Step 1: Replace all occurrences of x(n,m) with .{n,m}
        let mut regex_pattern = re_repeats.replace_all(pattern, |caps: &regex::Captures| {
            let n = &caps[1];
            let m = if &caps[2] != "" { &caps[2] } else { n }; // Handle single number case
            format!(".{{{},{}}}", n, m)
        }).to_string();

        // Step 2: Replace remaining 'x' with '.'
        regex_pattern = regex_pattern.replace("x", ".");

        // Step 3: Remove dashes ('-') as they are not necessary in regex
        regex_pattern = regex_pattern.replace("-", "");

        // Step 4: Ensure character sets like [IL] are preserved as-is
        regex_pattern
    }
}
impl SequenceFilter for HasSequenceMotif {
    fn filter(&self, sequence: &Sequence) -> bool { self.motif_regex.is_match(&sequence.to_string(0)) }
}