mod cif;

use pyo3::{pymodule, PyResult, Python};
use pyo3::types::PyModule;
use pyo3::prelude::*;
pub use cif::PyCifLoop;
pub use cif::PyCifData;

#[pymodule]
fn my_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCifLoop>()?;
    m.add_class::<PyCifData>()?;
    Ok(())
}

// mod cif;

// use pyo3::prelude::*;
// use pyo3::types::PyModule;

// pub use cif::PyCifLoop;
// pub use cif::PyCifData;

// /// This function is called to initialize the Python module.
// #[pymodule]
// fn my_module(py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_class::<PyCifLoop>()?;
//     m.add_class::<PyCifData>()?;
//     Ok(())
// }
