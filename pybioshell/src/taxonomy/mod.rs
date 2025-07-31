use pyo3::prelude::*;

mod py_rank;
use py_rank::{PyRank};

mod py_taxonomy;
use py_taxonomy::{PyNode, PyTaxonomy};

pub fn init_submodule(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent_module.py();
    let m = PyModule::new(py, "taxonomy")?;
    m.add_class::<PyRank>()?;
    m.add_class::<PyNode>()?;
    m.add_class::<PyTaxonomy>()?;
    parent_module.add_submodule(&m);

    // 🔽 Add to sys.modules["bioshell.taxonomy"]
    let sys_modules = py.import("sys")?.getattr("modules")?;
    sys_modules.set_item("bioshell.taxonomy", m)?;

    Ok(())
}

