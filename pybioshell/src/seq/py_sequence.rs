use pyo3::prelude::*;
use pyo3::exceptions::PyIndexError;

use bioshell_seq::sequence::Sequence;


#[pyclass(name = "Sequence")]
pub struct PySequence {
    pub(crate) inner: Sequence,
}

#[pymethods]
impl PySequence {
    /// Create a new `Sequence` instance from a description and a sequence string.
    ///
    /// Args:
    ///     description (str): The FASTA header or description line.
    ///     sequence (str): The amino acid or nucleotide sequence.
    ///
    /// Example:
    ///     >>> seq = Sequence("myseq", "ACDEFGHIKLMNPQRSTVWY")
    #[new]
    pub fn new(description: &str, seq: &str) -> Self {
        PySequence { inner: Sequence::new(description, seq) }
    }

    /// The description (header line) of this sequence.
    #[getter]
    pub fn description(&self) -> &str {
        self.inner.description()
    }

    /// The first n characters of the description line of this [`Sequence`].
    pub fn description_n(&self, n:usize) -> &str {
        self.inner.description_n(n)
    }

    /// Return the first identifier of a sequence, extracted from the description line.
    ///
    /// If `if_sort` is true, the identifiers will be sorted by arbitrary importance,
    /// so the method always reports the same type of identifier (e.g., GenBank or RefSeq) regardless their order.
    pub fn first_id(&self, if_sort: bool) -> String { self.inner.first_id(if_sort) }

    /// Return all identifiers of a sequence formatted into a single string.
    ///
    /// If `if_sort` is true, the identifiers will be sorted by arbitrary importance.
    pub fn full_id(&self, if_sort: bool) -> String { self.inner.full_id(if_sort) }

    /// The sequence itself as a string of characters.
    #[getter]
    pub fn seq(&self) -> String {
        String::from_utf8_lossy(self.inner.seq()).to_string()
    }

    /// Return the length of the sequence (number of residues).
    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Return the character (residue) at the given position.
    ///
    /// Raises:
    ///     IndexError: If the position is out of range.
    pub fn char(&self, pos: usize) -> PyResult<char> {
        if pos >= self.inner.len() {
            Err(PyIndexError::new_err("position out of range"))
        } else {
            Ok(self.inner.char(pos))
        }
    }

    /// Return the full sequence as a formatted string with optional line width.
    ///
    /// Args:
    ///     line_width (int): Maximum number of characters per line. Set to 0 to disable wrapping.
    pub fn format(&self, line_width: usize) -> String {
        self.inner.to_string(line_width)
    }

    /// Return the string representation of this sequence (unwrapped).
    fn __str__(&self) -> String {
        self.inner.to_string(0)
    }
}