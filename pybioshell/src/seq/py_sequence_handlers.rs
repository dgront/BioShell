use pyo3::prelude::*;

use crate::seq::{PySequence, PyAlignment};

/// Formats a pairwise alignment in the `pairwise` format.
///
/// Given the query and template sequences, an alignment path, and formatting parameters, this function returns a string representation of the alignment in a human-readable format.
/// `seq_name_width` is the width of the sequence name column, and `alignment_width` is the maximum number of alignment columns per block.
#[pyfunction]
#[pyo3(signature = (query, target, alignment, seq_name_width = 12, alignment_width = 70))]
pub fn format_pairwise(query: &PySequence, target: &PySequence,
                       alignment: &PyAlignment, seq_name_width: usize, alignment_width: usize) -> String {

    bioshell_seq::alignment::format_pairwise(&query.inner, &target.inner, &alignment.path, seq_name_width, alignment_width)
}