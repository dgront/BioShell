use pyo3::prelude::*;

pub mod seq;

#[pymodule]
fn bioshell(m: &Bound<'_, PyModule>) -> PyResult<()> {
    seq::init_submodule(m)?;
    Ok(())
}
