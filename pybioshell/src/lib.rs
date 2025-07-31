use pyo3::prelude::*;

pub mod seq;
pub mod taxonomy;

#[pymodule]
fn bioshell(m: &Bound<'_, PyModule>) -> PyResult<()> {
    seq::init_submodule(m)?;
    taxonomy::init_submodule(m)?;
    Ok(())
}
