mod bioshell_cif;

use pyo3::{Python, PyResult, PyModule, Exceptions};
use pyo3::prelude::PyModule;

use crate::bioshell_cif::{};

pub use bioshell_cif::*;

#[pymodule]
#[pyo3(name = "cif")]
fn cif_module(py: Python, m: &PyModule) -> PyResult<()> {        
    m.add_class::<CifError>()?;
    m.add_class::<CifData>()?;
    m.add_class::<CifLoop>()?;
    m.add_class::<CifDataItem>()?;
    m.add_class::<CifDataValue>()?;
    m.add_class::<CifDataValue::CifDataValueString>()?;
    m.add_class::<CifDataValue::CifDataValueFloat>()?;
    m.add_class::<CifDataValue::CifDataValueInt>()?;
    m.add_class::<CifDataValue::CifDataValueList>()?;
    m.add_class::<CifDataValue::CifDataValueLoop>()?;
    m.add_class::<CifDataValue::CifDataValueBlock>()?;    
    Ok(())
}