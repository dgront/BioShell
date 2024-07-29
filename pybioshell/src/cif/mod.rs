use py03::{Python, PyResult, PyModule, Exceptions};
use pyo3::prelude::*;

use bioshell_cif::*;


#[pyclass(name="CifLoop")]
#[derive(Debug,Clone,Default)]
pub struct PyCifLoop(CifLoop);

#[pymethods]
impl PyCifLoop{
    #[new]
    pub fn new() -> Self {
        PyCifLoop(CifLoop::new())
    }

    fn add_column(&mut self, column: String){
        if column.is_empty(){pybioshell
    }

    fn column_names(&self) -> Vec<String>{
        self.0.columns().clone()
    }
    
    fn columns(&self) -> Vec<String>{
        self.0.columns().clone()
    }

    fn rows(&self) -> Vec<Vec<String>>{
        self.0.rows().clone()
    }

    fn __str__(&self) -> String{
        self.0.to_string()
    }

    fn columm_mame_contains(&self, column: String) -> bool{
        self.0.columns().contains(&column)
    }

}

#[pyclass(name="CifDataValue")]
#[derive(Debug,Clone)]
pub struct PyCifDataValue{
    value: CifDataValue
}
   

#[pymethods]
impl PyCifDataValue{
    #[new]
    fn new_string(value: String) -> Self{
        PyCifDataValue::CifDataValueString(value)
    }

    #[new]
    fn new_float(value: f64) -> Self{
        PyCifDataValue::CifDataValueFloat(value)
    }

    #[new]
    fn new_int(value: i64) -> Self{
        PyCifDataValue::CifDataValueInt(value)
    }

    #[new]
    fn new_list(value: Vec<String>) -> Self{
        PyCifDataValue::CifDataValueList(value)
    }

    #[new]
    fn new_loop(value: PyCifLoop) -> Self{
        PyCifDataValue::CifDataValueLoop(value)
    }

    #[new]
    fn new_block(value: PyCifData) -> Self{
        PyCifDataValue::CifDataValueBlock(value)
    }

    fn __str__(&self) -> String{
        match self{
            PyCifDataValue::CifDataValueString(value) => value.clone(),
            PyCifDataValue::CifDataValueFloat(value) => value.to_string(),
            PyCifDataValue::CifDataValueInt(value) => value.to_string(),
            PyCifDataValue::CifDataValueList(value) => format!("{:?}", value),
            PyCifDataValue::CifDataValueLoop(value) => value.to_string(),
            PyCifDataValue::CifDataValueBlock(value) => value.to_string(),
        }
    }
}

impl From<PyCifDataValue> for CifDataValue{
    fn from(py_cif_data_value: PyCifDataValue) -> Self{
        match py_cif_data_value{
            PyCifDataValue::CifDataValueString(value) => CifDataValue::String(value),
            PyCifDataValue::CifDataValueFloat(value) => CifDataValue::Float(value),
            PyCifDataValue::CifDataValueInt(value) => CifDataValue::Int(value),
            PyCifDataValue::CifDataValueList(value) => CifDataValue::List(value),
            PyCifDataValue::CifDataValueLoop(value) => CifDataValue::Loop(value.into()),
            PyCifDataValue::CifDataValueBlock(value) => CifDataValue::Block(value.0),
        }
    }
}

impl From<CifDataValue> for PyCifDataValue{
    fn from(cif_data_value: CifDataValue) -> Self{
        match cif_data_value{
            CifDataValue::String(value) => PyCifDataValue::CifDataValueString(value),
            CifDataValue::Float(value) => PyCifDataValue::CifDataValueFloat(value),
            CifDataValue::Int(value) => PyCifDataValue::CifDataValueInt(value),
            CifDataValue::List(value) => PyCifDataValue::CifDataValueList(value),
            CifDataValue::Loop(value) => PyCifDataValue::CifDataValueLoop(PyCifLoop::from(value)),
            CifDataValue::Block(value) => PyCifDataValue::CifDataValueBlock(PyCifData(value)),
        }
    }
}


#[pyclass(name="CifDataItem")]
#[derive(Debug,Clone)]
pub struct PyCifDataItem(CifDataItem);




#[pyclass(name="CifData")]
#[derive(Debug,Clone)]
pub struct PyCifData(CifData);

#[pymethods]
impl PyCifData{
    #[new]
    fn new() -> Self{
        PyCifData(CifData::new())
    }

    fn add_item(&mut self, item: String, value: CifDataValue){
        self.0.add_item(item, value);
    }

    fn add_loop(&mut self, cif_loop: PyCifLoop){
        self.0.add_loop(cif_loop.0);
    }

    fn items(&self) -> Vec<(String, CifDataValue)>{
        self.0.items().clone()
    }

    fn loops(&self) -> Vec<PyCifLoop>{
        self.0.loops().iter().map(|x| PyCifLoop(x.clone())).collect()
    }

    fn __str__(&self) -> String{
        self.0.to_string()
    }
}

