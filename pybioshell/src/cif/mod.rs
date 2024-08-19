use py03::{Python, PyResult, PyModule, Exceptions};
use pyo3::prelude::*;

use bioshell_cif::*;


#[pyclass(name="CifLoop")]
pub struct PyCifLoop{
    column_mames: Vec<String>,
    data_rows: Vec<Vec<String>>,
    previoous_row_incomplete: bool,
};

#[pymethods]
impl PyCifLoop{
    #[new]
    pub fn new() -> Self {
        let cols : Vec<String> = data_items_names.iter().map(|x| x.to_string()).collect();
       return PyCifLoop{columm_mames: cols, data_rows: Vec::new(), previoous_row_incomplete: false};
    }

    fn add_column(&mut self, column: String){
       if self.data_rows.len() > 0{
           self.previoous_row_incomplete = true;
           panic!("Cannot add column after adding data rows");
       }
       self.columm_mames.push(column_names.to_string());
    }

    fn add_data_row(&mut self, row: Vec<String>){
        if row.len() != self.columm_mames.len(){
            panic!("Row length does not match column length");
        }
        if self.previoous_row_incomplete{
            panic!("Cannot add data row after adding incomplete data row");
        }
        self.data_rows.push(row);
    }    
    
    fn column_names(&self)  -> impl Iterator<Item = &String> { 
        return self.columm_mames.iter();
    }     

    fn data_rows(&self) -> impl Iterator<Item = &Vec<String>> {
        return self.data_rows.iter();
    }
    pub fn count_rows(&self) -> usize { self.data_rows.len() }
 
    pub fn count_columns(&self) -> usize { self.column_names.len() }
   
    pub fn column_index(&self, data_name: &str) -> Option<usize> {
        self.column_names.iter().position(|r| r == data_name)
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

    pub fn add_loop(&mut self, loop_: PyCifLoop){
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.add_loop(loop_),
            _ => panic!("Cannot add loop to non block value")
        }
    }

    pub fn add_item(&mut self, item: String, value: CifDataValue){
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.add_item(item, value),
            _ => panic!("Cannot add item to non block value")
        }
    }

    pub fn data_items(&self) -> Vec<(String, CifDataValue)>{
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.items(),
            _ => panic!("Cannot get data items from non block value")
        }
    }

    pub fn data_items_mut(&mut self) -> Vec<(String, CifDataValue)>{
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.items_mut(),
            _ => panic!("Cannot get data items from non block value")
        }
    }
    
    pub fn loop_blocks(&mut self) -> impl DoubleEndedIterator<Item = &mut CifLoop> + '_ {
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.loops(),
            _ => panic!("Cannot get loop blocks from non block value")
        }
    }

    pub fn first_loop(&self) -> Option<PyCifLoop>{
        match self{
            PyCifDataValue::CifDataValueBlock(value) => value.loops().first().cloned(),
            _ => panic!("Cannot get first loop from non block value")
        }
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

