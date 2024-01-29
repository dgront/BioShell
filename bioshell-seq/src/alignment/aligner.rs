use std::fmt::{Display, Formatter};
use std::slice::Iter;
use crate::sequence::Sequence;

/// Represents possible moves on an alignment matrix
pub enum AlignmentStep {
    /// Gap symbol inserted into a query (the first) sequence, denoted as ``'-'``
    ///
    /// Horizontal move eats a single position of a template sequence and inserts a gap in a query
    Horizontal,
    /// Gap symbol introduced in a template (the second) sequence, denoted as ``'|'``
    ///
    /// Vertical move eats a single position of a query sequence and inserts a gap in a template
    Vertical,
    /// Match between a query and a template symbols, denoted as ``'*'``
    Match
}

impl Display for AlignmentStep {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            AlignmentStep::Horizontal => { write!(f, "{}", "-")? }
            AlignmentStep::Vertical => { write!(f, "|")? }
            AlignmentStep::Match => { write!(f, "*")? }
        }
        Ok(())
    }
}

pub struct AlignmentPath {
    path: Vec<AlignmentStep>
}

impl AlignmentPath {
    pub fn from_attrs(path: Vec<AlignmentStep>) -> AlignmentPath { AlignmentPath{path} }

    pub fn iter(&self) -> Iter<'_, AlignmentStep> { self.path.iter() }
}


impl Display for AlignmentPath {
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

pub fn aligned_strings(alignment: &AlignmentPath, query: &str, template: &str, gap_symbol: char) -> (String, String) {
    let query_u8 = query.as_bytes();
    let template_u8 = template.as_bytes();
    let (ali_q, ali_t) = aligned_symbols(alignment, query_u8, template_u8, gap_symbol as u8);

    return (String::from_utf8(ali_q).unwrap(), String::from_utf8(ali_t).unwrap());
}

pub fn aligned_sequences(alignment: &AlignmentPath, query: &Sequence, template: &Sequence, gap_symbol: char) -> (Sequence, Sequence) {
    let (ali_q, ali_t) = aligned_symbols(alignment, query.as_u8(), template.as_u8(), gap_symbol as u8);

    (Sequence::from_attrs(query.description().to_string(), ali_q),
        Sequence::from_attrs(template.description().to_string(), ali_t))
}

pub trait Aligner<S> {
    fn align(&mut self) -> S;

}