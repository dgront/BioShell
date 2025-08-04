use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::basic::CompareOp;

use bioshell_taxonomy::Rank;

#[pyclass(name = "Rank")]
#[derive(Clone, Copy)]
pub struct PyRank { pub(crate) inner: Rank, }

#[pymethods]
impl PyRank {
    #[staticmethod]
    pub fn from_str(s: &str) -> PyResult<Self> {
        Ok(Self {
            inner: Rank::from_str(s),
        })
    }

    pub fn value(&self) -> u8 {
        self.inner as u8
    }

    pub fn name(&self) -> String {
        self.inner.to_string()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Rank::{}", self.inner.to_string())
    }

    fn __hash__(&self) -> u64 {
        self.inner as u64
    }

    fn __richcmp__(&self, other: PyRef<PyRank>, op: CompareOp) -> Py<PyAny> {
        Python::with_gil(|py| {
            let result = match op {
                CompareOp::Eq => self.inner == other.inner,
                CompareOp::Ne => self.inner != other.inner,
                CompareOp::Lt => (self.inner as u8) < (other.inner as u8),
                CompareOp::Le => (self.inner as u8) <= (other.inner as u8),
                CompareOp::Gt => (self.inner as u8) > (other.inner as u8),
                CompareOp::Ge => (self.inner as u8) >= (other.inner as u8),
            };
            result.into_py(py)
        })
    }
}

