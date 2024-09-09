use pyo3::prelude::*;
use pyo3::types::PyString;
use std::collections::HashMap;
use pyo3::{Python, PyClass, PyResult, pyclass, pymethods};



use bioshell_cif::CifLoop;
use bioshell_cif::CifData;


#[pyclass(name="CifLoop")]
pub struct PyCifLoop {
    data : CifLoop
}

#[pymethods]
impl PyCifLoop {
    #[new]
    pub fn new(data_item_names: Vec<String>) -> PyCifLoop {
        let cols: Vec<&str> = data_item_names.iter().map(|s| s.as_str()).collect();
        PyCifLoop {
            data: CifLoop::new(&cols),
        }
    }

    /// Add a new column
    pub fn add_column(&mut self, column_name: String) -> PyResult<()> {
        self.data.add_column(&column_name);
        Ok(())
    }

    /// Add a new row of data
    pub fn add_data_row(&mut self, row: Vec<String>) -> PyResult<()> {
        self.data.add_data_row(row);
        Ok(())
    }

    /// Get the list of rows
    pub fn rows(&self) -> Vec<Vec<String>> {
        self.data.rows().cloned().collect()
    }

    /// Get the list of column names
    pub fn column_names(&self) -> Vec<String> {
        self.data.column_names().cloned().collect()
    }

    /// Get the number of rows
    pub fn count_rows(&self) -> usize {
        self.data.count_rows()
    }

    /// Get the number of columns
    pub fn count_columns(&self) -> usize {
        self.data.count_columns()
    }

    /// Find the index of a column by name
    pub fn column_index(&self, data_name: String) -> Option<usize> {
        self.data.column_index(&data_name)
    }

    /// Check if a column name contains a specific substring
    pub fn column_name_contains(&self, substring: String) -> bool {
        self.data.column_name_contains(&substring)
    }

    /// Get a mutable reference to an entry (row and column) in the data
    pub fn entry_mut(&mut self, row_index: usize, data_name: String) -> Option<String> {
        self.data
            .entry_mut(row_index, &data_name)
            .map(|e| e.clone())  // Return a cloned version since we can't directly expose mutable references to Python
    }
}

#[pyclass(name = "CifData")]
pub struct PyCifData {
    data: CifData
}

#[pymethods]
impl PyCifData {
    #[new]
    pub fn new(name: &str) -> Self {
        PyCifData {
            data: CifData::new(name),
        }
    }

    pub fn name(&self) -> &str {
        &self.data.name()
    }

    pub fn add_item(&mut self, key: &str, value: &str) -> PyResult<()> {
        self.data.add_item(key, value);
        Ok(())
    }

    pub fn get_item(&self, key: &str) -> Option<String> {
        self.data.get_item(key)
    }

    pub fn data_items(&self) -> HashMap<String, String> {
        self.data.data_items().clone()
    }

    pub fn data_items_mut(&mut self) -> HashMap<String, String> {
        self.data.data_items_mut().clone()
    }

    //     pub fn add_loop(&mut self, a_loop: PyCifLoop) -> PyResult<()> {
    //     self.data.add_loop(a_loop.into());
    //     Ok(())
    // }

    // pub fn loop_blocks(&self) -> Vec<PyCifLoop> {
    //     self.data.loop_blocks().map(|l| PyCifLoop { loop_: l.clone() }).collect()
    // }

    // pub fn first_loop(&self, column: &str) -> Option<PyCifLoop> {
    //     self.data.first_loop(column).map(|l| PyCifLoop { loop_: l.clone() })
    // }
}