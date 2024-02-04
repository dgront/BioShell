use std::fmt::{Display, Formatter};
use std::slice::Iter;
use crate::sequence::Sequence;

/// Represents possible moves on an alignment matrix.
#[derive(Debug, PartialEq, Eq)]
pub enum AlignmentStep {
    /// Gap symbol inserted into a query (the first) sequence
    ///
    /// Horizontal step made on an alignment matrix eats a single position of a template sequence and inserts a gap in a query.
    /// It corresponds to the following alignment:
    /// ```text
    /// -
    /// A
    /// ```
    Horizontal,
    /// Gap symbol introduced in a template (the second) sequence
    ///
    /// Vertical step made on an alignment matrix eats a single position of a query sequence and inserts a gap in a template
    /// It corresponds to the following alignment:
    /// ```text
    /// A
    /// -
    /// ```
    Vertical,
    /// Match between a query and a template symbols corresponds to the following alignment:
    /// ```text
    /// I
    /// L
    /// ```
    Match
}

impl Display for AlignmentStep {
    /// Displays an alignment step as a single character
    ///
    /// The [`Horizontal`](AlignmentStep::Horizontal), [`Vertical`](AlignmentStep::Vertical)
    /// and [`Match`](AlignmentStep::Match) steps are displayed as ``'-'``, ``'|'`` and
    /// ``'*'``, respectively
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentStep::Horizontal => { write!(f, "{}", "-")? }
            AlignmentStep::Vertical => { write!(f, "|")? }
            AlignmentStep::Match => { write!(f, "*")? }
        }
        Ok(())
    }
}

impl TryFrom<u8> for AlignmentStep {
    type Error = &'static str;

    /// Tries to convert a `u8` value into an `AlignmentStep` variant.
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::alignment::AlignmentStep;
    /// assert_eq!(AlignmentStep::try_from(b'|').unwrap(), AlignmentStep::Vertical);
    /// assert_eq!(AlignmentStep::try_from(b'*').unwrap(), AlignmentStep::Match);
    /// ```
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'-' => Ok(AlignmentStep::Horizontal),
            b'|' => Ok(AlignmentStep::Vertical),
            b'*' => Ok(AlignmentStep::Match),
            _ => Err("Invalid value for AlignmentStep"),
        }
    }
}

/// Represents an abstract pairwise alignment.
///
/// An [AlignmentPath] object is an abstract definition of an alignment. It defines which position
/// in a query sequence should be aligned to a given position of a template and where gaps are located.
/// An alignment path is implemented a vector of [AlignmentStep]s taken on an alignment matrix to
/// align a query sequence with a template one.
pub struct AlignmentPath { path: Vec<AlignmentStep> }

impl AlignmentPath {
    /// Creates an [`AlignmentPath`] directly from steps
    pub fn from_attrs(path: Vec<AlignmentStep>) -> AlignmentPath { AlignmentPath{path} }

    /// Iterates over all steps of this path
    pub fn iter(&self) -> Iter<'_, AlignmentStep> { self.path.iter() }
}


impl TryFrom<&str> for AlignmentPath {
    type Error = &'static str;

    /// Tries to convert a `u8` value into an `AlignmentStep` variant.
    ///
    /// Each character of a given string is converted to an [AlignmentStep] variant with
    /// [`AlignmentStep::try_from(s: u8)`](AlignmentStep::try_from())
    /// # Example
    /// ```
    /// use bioshell_seq::alignment::AlignmentPath;
    /// let path = AlignmentPath::try_from("**-**").unwrap();
    /// assert_eq!(path.to_string(), "**-**");
    /// ```
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let path: Result<Vec<AlignmentStep>, _> = s.chars().map(|c| AlignmentStep::try_from(c as u8)).collect();
        path.map(|path| AlignmentPath { path })
    }
}

impl Display for AlignmentPath {
    /// Displays this [AlignmentPath] as a single line string
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for step in &self.path {
            write!(f, "{}", step)?;
        }
        Ok(())
    }
}

pub fn aligned_symbols<S: Copy>(alignment: &AlignmentPath, query: &[S], template: &[S], gap_symbol: S) -> (Vec<S>, Vec<S>) {
    let mut aliged_q: Vec<S> = vec![];
    let mut aliged_t: Vec<S> = vec![];
    let mut q_iter = query.iter();
    let mut t_iter = template.iter();
    for step in alignment.iter() {
        match step {
            AlignmentStep::Horizontal => {
                aliged_q.push(gap_symbol);
                aliged_t.push(t_iter.next().unwrap().clone());
            }
            AlignmentStep::Vertical => {
                aliged_q.push(q_iter.next().unwrap().clone());
                aliged_t.push(gap_symbol);
            }
            AlignmentStep::Match => {
                aliged_q.push(q_iter.next().unwrap().clone());
                aliged_t.push(t_iter.next().unwrap().clone());
            }
        }
    }
    return (aliged_q, aliged_t);
}

/// Expands two sequences into a sequence alignment based on an [AlignmentPath] object.
///
/// # Arguments
///
///  * `alignment` - [AlignmentPath] object is an abstract definition of an alignment
///  * `query` - query sequence (the first of the two aligned)
///  * `template` - template sequence (the second of the two aligned)
///  * `gap_symbol` - use ``'-'`` or  ``'_'``
///
/// The given [AlignmentPath] object defines an alignment, i.e. which letter of a query sequence should be aligned
/// to a given letter of a template and where gaps are located
///
/// # Example
/// ```
/// use bioshell_seq::alignment::{aligned_strings, AlignmentPath};
/// let alignment = AlignmentPath::try_from("**-**").unwrap();
/// let (aligned_query, aligned_template) = aligned_strings(&alignment, "ALIV", "ALRIV", '-');
/// assert_eq!(aligned_query, "AL-IV");
/// assert_eq!(aligned_template, "ALRIV");
/// ```
pub fn aligned_strings(alignment: &AlignmentPath, query: &str, template: &str, gap_symbol: char) -> (String, String) {
    let query_u8 = query.as_bytes();
    let template_u8 = template.as_bytes();
    let (ali_q, ali_t) = aligned_symbols(alignment, query_u8, template_u8, gap_symbol as u8);

    return (String::from_utf8(ali_q).unwrap(), String::from_utf8(ali_t).unwrap());
}

/// Expands two sequences into a sequence alignment based on an [AlignmentPath] object.
///
/// For instance, when the two (unaligned) sequences are ``"ALIV"`` and ``"ALRIV"`` and the alignment
/// path is ``"**-**"``, the expanded (aligned sequences) become ``"AL-IV"`` and ``"ALRIV"``
///
/// # Arguments
///
///  * `alignment` - [AlignmentPath] object is an abstract definition of an alignment
///  * `query` -  query sequence (the first of the two aligned)
///  * `template` - template sequence (the second of the two aligned)
///  * `gap_symbol` - use ``'-'`` or  ``'_'``
///
/// The given [AlignmentPath] object defines an alignment, i.e. which letter of a query sequence should be aligned
/// to a given letter of a template and where gaps are located
///
/// # Example
/// ```
/// use bioshell_seq::alignment::{aligned_sequences, AlignmentPath};
/// use bioshell_seq::sequence::Sequence;
/// let alignment = AlignmentPath::try_from("**-**").unwrap();
/// let query = Sequence::from_str("query","ALIV");
/// let template = Sequence::from_str("template","ALRIV");
/// let (aligned_query, aligned_template) = aligned_sequences(&alignment, &query, &template, '-');
/// assert_eq!(aligned_query.to_string(), "AL-IV");
/// assert_eq!(aligned_template.to_string(), "ALRIV");
/// ```
pub fn aligned_sequences(alignment: &AlignmentPath, query: &Sequence, template: &Sequence, gap_symbol: char) -> (Sequence, Sequence) {
    let (ali_q, ali_t) = aligned_symbols(alignment, query.as_u8(), template.as_u8(), gap_symbol as u8);

    (Sequence::from_attrs(query.description().to_string(), ali_q),
        Sequence::from_attrs(template.description().to_string(), ali_t))
}
