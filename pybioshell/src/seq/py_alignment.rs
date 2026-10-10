use pyo3::prelude::*;
use pyo3::exceptions::{PyValueError};

use bioshell_seq::alignment::{GlobalAligner, LocalAlignment, AlignmentPath, aligned_strings };
use bioshell_seq::scoring::SequenceSimilarityScore;

use crate::seq::PySequence;
use crate::seq::PySubstitutionMatrix;

//// Represents a sequence alignment.
#[pyclass(name = "Alignment")]
pub struct PyAlignment {
    score: i32,
    query_start: usize,
    template_start: usize,
    aligned_query: String,
    aligned_template: String,
}

impl PyAlignment {

    /// Creates a new `PyAlignment` instance.
    pub fn new(score: i32, aligned_query: String, aligned_template: String,
               query_start: usize, template_start: usize) -> Self {

        Self { score, query_start, template_start, aligned_query, aligned_template }
    }
}

#[pymethods]
impl PyAlignment {
    #[getter]
    fn score(&self) -> i32 { self.score }

    #[getter]
    fn query_start(&self) -> usize { self.query_start }

    #[getter]
    fn template_start(&self) -> usize { self.template_start }

    #[getter]
    fn aligned_query(&self) -> String { self.aligned_query.clone() }

    #[getter]
    fn aligned_template(&self) -> String { self.aligned_template.clone() }
}

#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn global_alignment(query: &PySequence, target: &PySequence,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {
    global_alignment_str(&query.seq(), &target.seq(), matrix, gap_open, gap_extend)
}

#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn global_alignment_str(query: &str, target: &str,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {

    let scoring = SequenceSimilarityScore::for_strings(query, target, matrix.matrix.clone());

    let mut aligner = GlobalAligner::new(query.len().max(target.len()));

    let score = aligner.align(&scoring, gap_open, gap_extend);
    let path = aligner.backtrace();
    let (aligned_query, aligned_template) = aligned_strings(&path, query, target, '-');
    Ok(PyAlignment::new(score, aligned_query, aligned_template, 0, 0))
}

#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open=-10, gap_extend=-1))]
pub fn local_alignment(query: &PySequence, target: &PySequence,
                        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {
    global_alignment_str(&query.seq(), &target.seq(), matrix, gap_open, gap_extend)
}

#[pyfunction]
#[pyo3(signature = (query, target, matrix, gap_open = - 10, gap_extend = - 1))]
pub fn local_alignment_str(query: &str, target: &str,
        matrix: &PySubstitutionMatrix, gap_open: i32, gap_extend: i32) -> PyResult<PyAlignment> {

    let scoring = SequenceSimilarityScore::for_strings(query, target, matrix.matrix.clone());
    let mut aligner = LocalAlignment::new(query.len().max(target.len()));
    let score = aligner.align(&scoring, gap_open, gap_extend);
    let (path, query_start, target_start) = aligner.backtrace();
    let (aligned_query, aligned_template) = aligned_strings(&path, query, target, '-');

    Ok(PyAlignment::new(score, aligned_query, aligned_template, query_start, target_start))
}