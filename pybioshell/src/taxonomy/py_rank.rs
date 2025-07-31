use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

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
}

