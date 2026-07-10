use std::fs::File;
use std::io::{BufReader};
use pyo3::prelude::*;
use pyo3::exceptions::{PyStopIteration, PyValueError};
use bioshell_seq::sequence::{FastaIterator, FastaParsingMode};
use crate::seq::PySequence;

#[doc = include_str!("../../docs/seq/fasta_iterator.rst")]
#[pyclass(name = "FastaIterator")]
pub struct PyFastaIterator {
    inner: FastaIterator<BufReader<File>>,
}

fn parsing_mode_from_str(parsing_mode: &str) -> Result<FastaParsingMode, PyErr> {

    match parsing_mode.trim().to_ascii_lowercase().as_str() {
        "raw" | "none" => Ok(FastaParsingMode::Raw),
        "protein" | "clean-protein" | "clean_protein" => Ok(FastaParsingMode::CleanProtein),
        "nucleic" | "clean-nucleic" | "clean_nucleic" => Ok(FastaParsingMode::CleanNucleic),
        "protein*" | "clean-protein*" | "clean_protein*"| "clean-protein-stop" | "clean_protein_stop" => Ok(FastaParsingMode::CleanProteinStop),

        _ => Err(PyValueError::new_err(format!(
            "invalid parsing_mode {parsing_mode:?}; expected one of: 'raw', 'protein', 'protein*', nucleic"
        ))),
    }
}


#[pymethods]
impl PyFastaIterator {
    /// Create a new iterator from a FASTA file `path`.
    ///
    /// Depending on the ``parsing_mode``, the file will be parsed:
    ///   - `"raw"` - as is, without any quality control
    ///   - `"protein"` - characters other than denoting amino acids will be removed
    ///   - `"protein*"` - as in  `"protein"` but stop codon denoted as `*` is allowed
    ///   - `"nucleic"` - characters other than denoting nucleic acids will be removed
    ///
    /// Returns:
    ///     FastaIterator: An iterator over sequences.
    #[new]
    pub fn new(path: &str, parsing_mode: &str) -> PyResult<Self> {
        let file = File::open(path).map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to open file: {}", e)))?;
        let reader = BufReader::new(file);
        Ok(PyFastaIterator {
            inner: FastaIterator::new(reader, parsing_mode_from_str(parsing_mode)?),
        })
    }

    /// Return the iterator itself.
    fn __iter__(slf: PyRefMut<'_, Self>) -> PyRefMut<'_, Self> {
        slf
    }

    /// Return the next sequence from the FASTA file.
    fn __next__(&mut self) -> PyResult<PySequence> {
        match self.inner.next() {
            Some(result) => result
                .map(|seq| PySequence { inner: seq })
                .map_err(|err| PyValueError::new_err(err.to_string())),
            None => Err(PyStopIteration::new_err("End of FASTA file")),
        }
    }
}
