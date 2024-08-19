// pub mod pybioshell;

// use pyo3::prelude::*;
// use pyo3::wrap_pyfunction;

// use crate::pybioshell::PyCifLoop;

#[pymodule]
fn pybioshell(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCifLoop>()?;
    Ok(())
}


