use pyo3::prelude::*;

mod py_sequence;
use py_sequence::{PySequence};

mod py_fasta_iterator;
use py_fasta_iterator::{PyFastaIterator};

pub fn init_submodule(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent_module.py();
    let m = PyModule::new(py, "seq")?;
    m.add_class::<PySequence>()?;
    m.add_class::<PyFastaIterator>()?;
    parent_module.add_submodule(&m);

    // 🔽 Add to sys.modules["bioshell.seq"]
    let sys_modules = py.import("sys")?.getattr("modules")?;
    sys_modules.set_item("bioshell.seq", m)?;

    Ok(())
}

