use py03::{Python, PyResult, PyModule, exceptions};
use pyo3::prelude::*;

use bioshell_cif::*;

#[pyclass(name="CifLoop")]
#[derive(Debug,clone)]
pub struct PyCifLoop(CifLoop);


#[pymethods]
impl PyCifLoop{
    #[new]
    fn new() -> Self{
        PyCifLoop(CifLoop::new())
    }

    fn add_column(&mut self, column: String){
        self.0.add_column(column);
    }

    fn add_row(&mut self, row: Vec<String>){
        self.0.add_row(row);
    }

    fn columns(&self) -> Vec<String>{
        self.0.columns().clone()
    }

    fn rows(&self) -> Vec<Vec<String>>{
        self.0.rows().clone()
    }
}

