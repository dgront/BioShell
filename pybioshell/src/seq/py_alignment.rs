use pyo3::prelude::*;
use bioshell_seq::alignment::{GlobalAligner, LocalAlignment, AlignmentPath, aligned_strings, aligned_sequences};
use bioshell_seq::scoring::SequenceSimilarityScore;

use crate::seq::PySequence;
use crate::seq::PySubstitutionMatrix;

//// Represents an abstract pairwise sequence alignment.
#[pyclass(name = "Alignment")]
pub struct PyAlignment {
    score: i32,
    query_start: usize,
    template_start: usize,
    path: AlignmentPath,
}

impl PyAlignment {

    /// Creates a new `PyAlignment` instance.
    pub fn new(score: i32, query_start: usize, template_start: usize, path: AlignmentPath) -> Self {

        Self { score, query_start, template_start, path }
    }
}

#[pymethods]
impl PyAlignment {
    /// Score of the alignment, which is the sum of substitution scores and gap penalties.
    #[getter]
    fn score(&self) -> i32 { self.score }

    /// Starting position of the aligned query sequence - zero-based index.
    ///
    /// In the case of global alignment, this is always 0. In the case of local alignment, this is the index of the first aligned residue in the query sequence.
    #[getter]
    fn query_start(&self) -> usize { self.query_start }

    /// Starting position of the aligned template sequence - zero-based index.
    ///
    /// In the case of global alignment, this is always 0. In the case of local alignment, this is the index of the first aligned residue in the template sequence.
    #[getter]
    fn template_start(&self) -> usize { self.template_start }

    /// Returns the aligned query and template sequences as a tuple of Sequence objects.
    fn aligned_sequences(&self, full_query: &PySequence, full_template: &PySequence) -> (PySequence, PySequence) {

        let (ali_q, ali_t) = aligned_sequences(&self.path, &full_query.inner, &full_template.inner, '-');
        (PySequence { inner: ali_q }, PySequence { inner: ali_t })
    }

    /// Returns the aligned query and template sequences as a tuple of strings.
    fn aligned_str(&self, full_query: &str, full_template: &str) -> (String,String) {
        aligned_strings(&self.path, full_query, full_template, '-')
    }
}

/// Aligns two sequences globally using the Needleman-Wunsch algorithm with affine gap penalties.
///
/// Returns a `Alignment` object, which allows to retrieve the aligned sequences.
#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn global_alignment(query: &PySequence, target: &PySequence,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {
    global_alignment_str(&query.seq(), &target.seq(), matrix, gap_open, gap_extend)
}


/// Aligns two sequences globally using the Needleman-Wunsch algorithm with affine gap penalties.
///
/// Returns a `Alignment` object, which allows to retrieve the aligned sequences.
#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn global_alignment_str(query: &str, target: &str,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {

    let scoring = SequenceSimilarityScore::for_strings(query, target, matrix.matrix.clone());

    let mut aligner = GlobalAligner::new(query.len().max(target.len()));

    let score = aligner.align(&scoring, gap_open, gap_extend);
    let path = aligner.backtrace();
    Ok(PyAlignment::new(score, 0, 0, path))
}

/// Aligns two sequences locally using the Smith-Waterman algorithm with affine gap penalties.
///
/// Returns a `Alignment` object, which allows to retrieve the aligned sequences.
#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn local_alignment(query: &PySequence, target: &PySequence,
                        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {
    global_alignment_str(&query.seq(), &target.seq(), matrix, gap_open, gap_extend)
}

/// Aligns two sequences locally using the Smith-Waterman algorithm with affine gap penalties.
///
/// Returns a `Alignment` object, which allows to retrieve the aligned sequences.
#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open = - 10, gap_extend = - 1))]
pub fn local_alignment_str(query: &str, target: &str,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {

    let scoring = SequenceSimilarityScore::for_strings(query, target, matrix.matrix.clone());
    let mut aligner = LocalAlignment::new(query.len().max(target.len()));
    let score = aligner.align(&scoring, gap_open, gap_extend);
    let (path, query_start, target_start) = aligner.backtrace();

    Ok(PyAlignment::new(score, query_start, target_start, path))
}