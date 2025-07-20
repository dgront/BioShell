use pyo3::prelude::*;
use bioshell_seq::sequence::Sequence;

#[pyclass(name = "Sequence")]
pub struct PySequence {
    inner: Sequence,
}

#[pymethods]
impl PySequence {
    #[new]
    pub fn new(description: &str, seq: &str) -> Self {
        PySequence { inner: Sequence::new(description, seq) }
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

pub fn init_submodule(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent_module.py();
    let m = PyModule::new(py, "seq")?;
    m.add_class::<PySequence>()?;
    parent_module.add_submodule(&m);

    // 🔽 Add to sys.modules["bioshell.seq"]
    let sys_modules = py.import("sys")?.getattr("modules")?;
    sys_modules.set_item("bioshell.seq", m)?;

    Ok(())
}

