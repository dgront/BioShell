use pyo3::prelude::*;

mod py_sequence;
use py_sequence::{PySequence};

mod py_sequence_id;
use py_sequence_id::{PySeqId, parse_sequence_id};

mod py_fasta_iterator;
use py_fasta_iterator::{PyFastaIterator};
mod py_substitution_matrix;
use py_substitution_matrix::{PySubstitutionMatrix};
mod py_alignment;
use py_alignment::{PyAlignment, global_alignment, local_alignment, global_alignment_str, local_alignment_str};

mod py_sequence_handlers;
use py_sequence_handlers::{format_pairwise};



pub fn init_submodule(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent_module.py();
    let m = PyModule::new(py, "seq")?;
    m.add_class::<PySeqId>()?;
    m.add_class::<PySequence>()?;
    m.add_class::<PyFastaIterator>()?;
    m.add_class::<PySubstitutionMatrix>()?;
    m.add_class::<PyAlignment>()?;

    m.add_function(wrap_pyfunction!(parse_sequence_id, &m)?)?;
    m.add_function(wrap_pyfunction!(global_alignment, &m)?)?;
    m.add_function(wrap_pyfunction!(local_alignment, &m)?)?;
    m.add_function(wrap_pyfunction!(global_alignment_str, &m)?)?;
    m.add_function(wrap_pyfunction!(local_alignment_str, &m)?)?;
    m.add_function(wrap_pyfunction!(format_pairwise, &m)?)?;
    parent_module.add_submodule(&m)?;

    // 🔽 Add to sys.modules["bioshell.seq"]
    let sys_modules = py.import("sys")?.getattr("modules")?;
    sys_modules.set_item("bioshell.seq", m)?;

    Ok(())
}

