use crate::sequence::{count_residue_type, Sequence};

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

/// Always returns `true`; serves a the defafult filter that pass every sequence
pub struct AlwaysTrue;

impl SequenceFilter for AlwaysTrue {
    fn filter(&self, _: &Sequence) -> bool { true }
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