use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader};
use log::{debug, info};

pub struct CifData {
    name: String,
    entries: HashMap<String, String>,
    loops: Vec<CifLoop>
}

pub struct CifLoop {
    column_names: Vec<String>,
    data_rows: Vec<Vec<String>>
}

impl CifLoop {

    /// Add a new column to this loop block.
    ///
    /// Adding columns is possible only before any data is inserted; once any data has been inserted,
    /// this method will panic.
    pub fn add_column(&mut self, column_name: String) {
        if self.data_rows.len() > 0 { panic!("Addempted column insertion for a loop-block that already contains some data!"); }
        self.column_names.push(column_name);
    }

    /// Add a new row of data.
    ///
    /// The provided row of data must contain the same number of entries as the number of columns
    /// in this loop-block; otherwise this method will panic.
    pub fn add_data_row(&mut self, row: Vec<String>) {
        if self.column_names.len() != row.len() {
            // panic!("Provided row of data doesn't match the number of columns!\nOffending input was:{:?}", &row);
        }
        self.data_rows.push(row);
    }
}

impl Display for CifLoop {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "loop_").ok();
        for coln in &self.column_names {
            writeln!(f, "{}", coln).ok();
        }
        for row in &self.data_rows {
            for val in row {
                write!(f, " {}",val).ok();
            }
            writeln!(f, "").ok();
        }
        Ok(())
    }
}

impl CifData {
    pub fn new(name:String) -> CifData {
        return CifData{
            name, entries: HashMap::new(), loops: vec![]
        };
    }

    pub fn insert(&mut self, key: String, value: String) {
        self.entries.insert(key, value);
    }

    /// Add a new empty loop to this block
    pub fn new_loop(&mut self) { self.loops.push(CifLoop{ column_names: vec![], data_rows: vec![] })}

    /// Add a given loop to this block
    pub fn add_loop(&mut self, a_loop: CifLoop) { self.loops.push(a_loop)}

    /// Get an iterator of references to loop-blocks.
    pub fn loop_blocks(&self) -> impl DoubleEndedIterator<Item = &CifLoop> + '_ {
        self.loops.iter()
    }

    /// Get an iterator of mutable references to loop-blocks.
    pub fn loop_blocks_mut(&mut self) -> impl DoubleEndedIterator<Item = &mut CifLoop> + '_ {
        self.loops.iter_mut()
    }
}

impl Display for CifData {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (key, val) in &self.entries {
            writeln!(f, "{} {}", key, val).ok();
        }
        writeln!(f,"").ok();
        for a_loop in &self.loops {
            writeln!(f, "{}",a_loop).ok();
        }
        Ok(())
    }
}

pub fn read_cif_file(input_fname: &str) -> Result<Vec<CifData>, io::Error> {

    let file = match File::open(input_fname) {
        Ok(file) => file,
        Err(e) => return Err(e)
    };
    return Ok(read_cif_buffer(&mut BufReader::new(file)));
}

pub fn read_cif_buffer<R: BufRead>(buffer: &mut R) -> Vec<CifData> {

    let mut data_blocks: Vec<CifData> = vec![];
    let mut current_block: Option<CifData> = None;
    let mut current_loop = CifLoop{ column_names: vec![], data_rows: vec![] };
    let mut is_loop_open: bool = false;

    for line in buffer.lines() {
        if let Some(line_ok) = line.ok() {
            let ls = line_ok.trim();
            // ---------- skip empty content and comments
            if ls.len() == 0 || ls.starts_with("#") { continue; }

            if ls.starts_with("data_") {        // --- start a new data block
                let mut parts = ls.splitn(2, '_');
                let block_name = parts.nth(1).unwrap().to_string();
                info!("new data block found: {}",&block_name);
                data_blocks.push(CifData::new(block_name));
            } else if ls.starts_with('_') {
                if is_loop_open {         // --- Add a key to the current loop block
                    current_loop.add_column(ls.to_string());
                } else {
                    if data_blocks.len() == 0 { panic!("Found data entries outside any data block!")}
                    let mut key_val = ls.split_whitespace();
                    let key = key_val.nth(0).unwrap();
                    let mut val: String = "".to_string();
                    for s in &mut key_val { val += &s.to_string(); }
                    data_blocks.last_mut().unwrap().insert(key.to_string(), val);
                }
            } else if ls.starts_with("loop_") {
                if data_blocks.len() == 0 { panic!("Found data loop outside any data block!")}
                if is_loop_open {
                    data_blocks.last_mut().unwrap().add_loop(current_loop);
                }
                current_loop = CifLoop{ column_names: vec![], data_rows: vec![] };
                is_loop_open = true;
            } else {
                current_loop.add_data_row(ls.split_whitespace().map(|s| s.to_string()).collect());
            }
        }
    }

    return data_blocks;
}