mod cif;
mod seq; // For pybioshell/src/seq/mod.rs

use pyo3::{pymodule, PyResult, Python};
use pyo3::types::PyModule;
use pyo3::prelude::*;

use pyo3::wrap_pyfunction; // For pybioshell/src/seq/mod.rs
pub use cif::PyCifLoop; 
pub use cif::PyCifData;

pub use seq::MonomerType; // For pybioshell/src/seq/mod.rs
pub use seq::get_monomer_type_value; // For pybioshell/src/seq/mod.rs

#[pymodule]
fn my_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyCifLoop>()?;
    m.add_class::<PyCifData>()?;
    Ok(())
}

#[pymodule] // For pybioshell/src/seq/mod.rs
fn residue_types(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MonomerType>()?;
    m.add_function(wrap_pyfunction!(get_monomer_type_value, m)?)?;
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
