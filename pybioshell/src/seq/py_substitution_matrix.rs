use std::path::Path;
use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};

use pyo3::exceptions::{PyValueError, PyIOError};
use pyo3::prelude::*;

#[pyclass(name = "SubstitutionMatrix")]
/// Represents a substitution matrix for scoring amino acid similarity.
/// The matrix can be loaded from a predefined list of matrices (e.g., BLOSUM62, PAM250)
/// or from an external file in the NCBI format.
pub struct PySubstitutionMatrix {
    name: String, //  the matrix description
    pub(crate) matrix: SubstitutionMatrix
}

fn matrix_kind(name: &str) -> PyResult<SubstitutionMatrixList> {
    match name.to_ascii_uppercase().as_str() {
        "BLOSUM45" => Ok(SubstitutionMatrixList::BLOSUM45),
        "BLOSUM62" => Ok(SubstitutionMatrixList::BLOSUM62),
        "BLOSUM80" => Ok(SubstitutionMatrixList::BLOSUM80),
        "PAM30" => Ok(SubstitutionMatrixList::PAM30),
        "PAM70" => Ok(SubstitutionMatrixList::PAM70),
        "PAM120" => Ok(SubstitutionMatrixList::PAM120),
        "PAM250" => Ok(SubstitutionMatrixList::PAM250),
        _ => Err(PyValueError::new_err(
            format!("Unknown substitution matrix: {name}")
        )),
    }
}


#[pymethods]
impl PySubstitutionMatrix {
    #[new]
    #[pyo3(signature = (name="BLOSUM62"))]
    fn new(name: &str) -> PyResult<Self> {
        let kind = matrix_kind(name)?;
        Ok(Self {
            name: name.to_ascii_uppercase(),
            matrix: SubstitutionMatrix::load(kind),
        })
    }

    /// Load a substitution matrix from an NCBI-formatted file.
    #[staticmethod]
    fn from_file(file_name: &str) -> PyResult<Self> {
        let matrix = SubstitutionMatrix::ncbi_matrix_from_file(file_name)
            .map_err(|e| PyIOError::new_err(e.to_string()))?;

        let name = Path::new(file_name)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(file_name)
            .to_string();

        Ok(Self { name, matrix })
    }

    #[getter]
    /// Returns the name of the substitution matrix (e.g., "BLOSUM62", "PAM250").
    fn name(&self) -> &str { &self.name }

    /// Returns the score for aligning two amino acids represented by their one-letter codes.
    fn score(&self, a: char, b: char) -> PyResult<i32> {

        Ok(self.matrix.score_by_aa(a as u8, b as u8))
    }

    fn __repr__(&self) -> String {
        format!("SubstitutionMatrix('{}')", self.name)
    }
}

