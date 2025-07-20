use std::fs::File;
use std::io::{BufReader};
use pyo3::prelude::*;
use pyo3::exceptions::PyStopIteration;
use bioshell_seq::sequence::FastaIterator;
use crate::seq::PySequence;

/// An iterator over sequences in a FASTA file.
///
/// Yields `Sequence` objects one by one.
#[pyclass(name = "FastaIterator")]
pub struct PyFastaIterator {
    inner: FastaIterator<BufReader<File>>,
}

#[pymethods]
impl PyFastaIterator {
    /// Create a new iterator from a FASTA file path.
    ///
    /// Args:
    ///     path (str): Path to the FASTA file.
    ///
    /// Returns:
    ///     PyFastaIterator: An iterator over sequences.
    #[new]
    pub fn new(path: &str) -> PyResult<Self> {
        let file = File::open(path).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open file: {}", e)))?;
        let reader = BufReader::new(file);
        Ok(PyFastaIterator {
            inner: FastaIterator::new(reader),
        })
    }

    /// Return the iterator itself.
    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    /// Return the next sequence from the FASTA file.
    ///
    /// Raises:
    ///     StopIteration: When there are no more sequences.
    fn __next__(&mut self) -> PyResult<PySequence> {
        match self.inner.next() {
            Some(seq) => Ok(PySequence { inner: seq }),
            None => Err(PyStopIteration::new_err("End of FASTA file")),
        }
    }
}
