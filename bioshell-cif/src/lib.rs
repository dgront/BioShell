//! Reads and writes data in CIF format.
//!
//! A CIF file stores **data blocks**. Each data block in turn contains name-value data items and loop blocks.
//! A single file may contain more than one data block. In addition, the CIF format defines **loop blocks**,
//! which store tabulated data.
//!
//! # Example CIF-formatted file:
//! ``` text
//! data_some_name
//! _name_1            value_1
//! _name_2            value_2
//!
//! loop_
//! _first_column
//! _second_column
//! 'value A' 1
//! 'value B' 2
//! 'value C' 2
//! ```
//!
//! This example CIF entry contains a single block, named ``some_name`` (the mandatory ``data_``
//! prefix is not considered a part of that name). That block, loaded as a [CifData]
//! struct, holds two key-value entries:
//! ``_name_1:value_1`` and ``_name_2:value_2``,  followed by a loop block.
//! Data items stored by a loop block are loaded into a [CifLoop] struct.
//!
//! The official specification of the CIF format can be found on
//! [this page](https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax)
//!
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Lines};
use std::time::Instant;
use log::{debug, info};

use bioshell_io::split_into_strings;


/// Returns true if a given file is in CIF format.
///
/// This function simply tests whether the first data line of a given file starts with ``data_``.
/// Otherwise it returns ``false``. When the file can't be open returns I/O error..
pub fn is_cif_file(file_path: &str) -> io::Result<bool> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        if !line.trim().is_empty() && !line.starts_with('#') {
            return Ok(line.trim().starts_with("data_"));
        }
    }

    return Ok(false);
}

/// Represents a single `data_` block of a CIF file.
///
/// A single data block may contain entries given as key-value pairs as well as loop blocks.
pub struct CifData {
    name: String,
    data_items: HashMap<String, String>,
    loops: Vec<CifLoop>
}

/// Represents a single `_loop` of a CIF file.
///
/// # Example
///
/// The following example shows, how to build a [CifLoop] using its API:
/// ```
/// use bioshell_cif::CifLoop;
/// // --- create an empty data loop with four columns
/// let mut data_loop = CifLoop::new(&["_atom_site_label", "_atom_site_Cartn_x",
///         "_atom_site_Cartn_y", "_atom_site_Cartn_z"]);
/// // --- append two rows of values
/// data_loop.add_data_row(vec!["O1", "4.154", "5.699", "3.026"].iter().map(|&s| s.to_string()).collect());
/// data_loop.add_data_row(vec!["C2", "5.630", "5.087", "3.246"].iter().map(|&s| s.to_string()).collect());
/// // --- modify one of the entries of the table
/// *data_loop.entry_mut(1, "_atom_site_Cartn_z").unwrap() = "4.246".to_string();
/// // --- print the data loop block
/// println!("{}", data_loop);
/// # let out = format!("{}", data_loop);
/// # let exp = "loop_
/// # _atom_site_label
/// # _atom_site_Cartn_x
/// # _atom_site_Cartn_y
/// # _atom_site_Cartn_z
/// # O1 4.154 5.699 3.026
/// # C2 5.630 5.087 4.246
///
/// # ";
/// # assert_eq!(out, exp);
/// ```
/// The script should print the following result:
/// ```text
/// loop_
/// _atom_site_label
/// _atom_site_Cartn_x
/// _atom_site_Cartn_y
/// _atom_site_Cartn_z
/// O1 4.154 5.699 3.026
/// C2 5.630 5.087 4.246
/// ```
pub struct CifLoop {
    column_names: Vec<String>,
    data_rows: Vec<Vec<String>>,
    previous_row_incomplete: bool
}

impl CifLoop {

    /// Creates an empty loop block with given columns.
    ///
    /// The newly created struct basically represents a table with named columns but with no data rows
    pub fn new(data_item_names: &[&str]) -> CifLoop {
        let cols: Vec<_> = data_item_names.iter().map(|e| e.to_string()).collect();
        return CifLoop{ column_names: cols, data_rows: vec![], previous_row_incomplete: false };
    }

    /// Add a new column to this loop block.
    ///
    /// Adding columns is possible only before any data is inserted; once any data has been inserted,
    /// this method will panic.
    pub fn add_column(&mut self, column_name: &str) {
        if self.data_rows.len() > 0 {
            panic!("Attempted column insertion for a loop-block that already contains some data!");
        }
        self.column_names.push(column_name.to_string());
    }

    /// Add a new row of data.
    ///
    /// The provided row of data must contain the same number of entries as the number of columns
    /// in this loop-block; otherwise this method will panic.
    pub fn add_data_row(&mut self, row: Vec<String>) {

        let n_columns = self.column_names.len();
        if self.column_names.len() != row.len() {
            if self.previous_row_incomplete {
                let last_vec: &mut Vec<String> = self.data_rows.last_mut().unwrap();
                last_vec.extend(row);
                if last_vec.len() == n_columns { self.previous_row_incomplete = false; }
                else if last_vec.len() > n_columns {
                    panic!("Provided row of data doesn't match the number of columns!\nThe previous row was incomplete; the combined rows are: {:?}", &last_vec);
                }
            } else {
                self.data_rows.push(row);
                self.previous_row_incomplete = true;
            }
        }
    }

    /// Non-mutable iterator over rows of this loop block.
    pub fn rows(&self) -> impl Iterator<Item = &Vec<String>> {
        return self.data_rows.iter();
    }

    /// Non-mutable iterator over names assigned to the columns of this loop.
    pub fn column_names(&self)  -> impl Iterator<Item = &String> { return self.column_names.iter(); }

    /// Counts rows of data stored by this loop
    pub fn count_rows(&self) -> usize { self.data_rows.len() }

    /// Counts columns (i.e. data items) stored by this loop
    pub fn count_columns(&self) -> usize { self.column_names.len() }

    /// Index of a column which holds values for a data item given its name
    pub fn column_index(&self, data_name: &str) -> Option<usize> {
        self.column_names.iter().position(|r| r == data_name)
    }

    /// Provides access to a data item from a given row of this loop
    ///
    /// This method allows change a single entry of this data loop
    pub fn entry_mut(&mut self, row_index: usize, data_name: &str) -> Option<&mut String> {
        return match self.column_index(data_name) {
            None => { None }
            Some(idx) => { Some(&mut self.data_rows[row_index][idx]) }
        };
    }
}

impl Display for CifLoop {
    /// Writes a [`CifLoop`](CifLoop) block in the CIF format.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_some_name
    /// loop_
    /// _first_column
    /// _second_column
    /// 'value A' 1
    /// 'value B' 2
    /// 'value C' 2
    /// #
    /// ";
    /// let mut reader = BufReader::new(cif_block.as_bytes());
    /// let data_blocks = read_cif_buffer(&mut reader);
    /// assert_eq!(data_blocks.len(), 1);
    /// assert_eq!(data_blocks[0].name(),"some_name");
    /// let a_loop = data_blocks[0].loop_blocks().next().unwrap();
    /// let out = format!("{}", a_loop);
    /// println!("{}", out);
    /// ```
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "loop_").ok();
        for coln in &self.column_names {
            writeln!(f, "{}", coln).ok();
        }
        for row in &self.data_rows {
            write!(f, "{}", &row[0]).ok();
            for val in &row[1..] {
                write!(f, " {}", val).ok();
            }
            writeln!(f, "").ok();
        }
        writeln!(f, "").ok();   // --- new line after a loop block
        Ok(())
    }
}

impl CifData {
    /// Creates a new CIF data structure with the given name.
    ///
    /// This method creates an empty struct with no data.
    ///
    /// # Example
    /// ```
    /// use bioshell_cif::CifData;
    /// let cif_block = CifData::new("blockname");
    /// assert_eq!(cif_block.name(), "blockname");
    /// ```
    pub fn new(name: &str) -> CifData {
        return CifData{
            name: name.to_string(), data_items: HashMap::new(), loops: vec![]
        };
    }

    /// name of this data block
    ///
    /// # Examples
    /// ```rust
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_first_block";
    /// let mut reader = BufReader::new(cif_block.as_bytes());
    /// let data_blocks = read_cif_buffer(&mut reader);
    /// assert_eq!(data_blocks.len(), 1);
    /// assert_eq!(data_blocks[0].name(),"first_block");
    /// ```
    pub fn name(&self) -> &str { &self.name }

    /// Add a given loop to this block
    pub fn add_loop(&mut self, a_loop: CifLoop) { self.loops.push(a_loop)}

    /// Read access to data items of this block
    pub fn data_items(&self) -> &HashMap<String, String> { &self.data_items }

    /// Mutable access to data items of this block
    pub fn data_items_mut(&mut self) -> &mut HashMap<String, String> { &mut self.data_items }

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
    /// Writes a [`CifData`](CifData) block in the CIF format.
    /// All loop-blocks contained in this block will also be displayed.
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let max_key_length = self.data_items.keys()
            .map(|key| key.len())
            .fold(0, |max_length, key_length| key_length.max(max_length));

        writeln!(f,"data_{}", self.name()).ok();
        for (key, val) in &self.data_items {
            writeln!(f, "{:width$} {}", key, val, width = max_key_length).ok();
        }
        writeln!(f,"").ok();
        for a_loop in &self.loops {
            writeln!(f, "{}",a_loop).ok();
        }
        Ok(())
    }
}

/// Reads a CIF-formatted file.
///
/// This function opens a file as a Returns a ``BufRead``, calls (``read_cif_buffer()``)[read_cif_buffer()]
/// and returs all the data blocks it found.
pub fn read_cif_file(input_fname: &str) -> Result<Vec<CifData>, io::Error> {

    info!("Loading a CIF file: {}", input_fname);

    let file = match File::open(input_fname) {
        Ok(file) => file,
        Err(e) => return Err(e)
    };
    return Ok(read_cif_buffer(&mut BufReader::new(file)));
}

/// Reads a CIF-formatted data from a buffer.
///
/// Returns a vector of all data blocks found.
pub fn read_cif_buffer<R: BufRead>(buffer: R) -> Vec<CifData> {

    let mut data_blocks: Vec<CifData> = vec![];
    let mut current_loop: Option<CifLoop> = None;
    // let mut current_loop = CifLoop{ column_names: vec![], data_rows: vec![] };
    let mut is_loop_open: bool = false;

    let start = Instant::now();
    let mut line_iter = buffer.lines();
    while let Some(line) = line_iter.next() {
        if let Some(line_ok) = line.ok() {
            let ls = line_ok.trim();
            // ---------- skip empty content and comments
            if ls.len() == 0 || ls.starts_with("#") {
                if is_loop_open {
                    data_blocks.last_mut().unwrap().add_loop(current_loop.unwrap());
                    current_loop = None;
                }
                is_loop_open = false;
                continue;
            }

            if ls.starts_with("data_") {        // --- start a new data block
                let mut parts = ls.splitn(2, '_');
                let block_name = parts.nth(1).unwrap();
                info!("new data block found: {}",&block_name);
                data_blocks.push(CifData::new(block_name));
            } else if ls.starts_with('_') {
                if is_loop_open {         // --- Add a key to the current loop block
                    if let Some(a_loop) = &mut current_loop {
                        a_loop.add_column(ls);
                    } else {
                        panic!("Attempt to add a column with no loop open");
                    }
                } else {
                    if data_blocks.len() == 0 { panic!("Found data entries outside any data block!")}
                    let last_block = data_blocks.last_mut().unwrap();
                    let key_val = split_into_strings(ls, false);
                    if key_val.len() != 2 {
                        let next_line = read_string("", &mut line_iter);
                        if ! is_quoted_string(&next_line.trim()) {
                            panic!("{}", format!("A single data item line should contain exactly two tokens: a key and its value; {} values found in the line: {}",
                                                 key_val.len(), &ls));
                        }
                        last_block.data_items_mut().insert(key_val[0].clone(), next_line.trim().to_string());
                    } else {
                        last_block.data_items_mut().insert(key_val[0].clone(), key_val[1].clone());
                    }
                }
            } else if ls.starts_with("loop_") {
                if data_blocks.len() == 0 { panic!("Found data loop outside any data block!")}
                if is_loop_open {
                    data_blocks.last_mut().unwrap().add_loop(current_loop.unwrap());
                }
                current_loop = Some(CifLoop{ column_names: vec![], data_rows: vec![], previous_row_incomplete: false });
                is_loop_open = true;
            } else {        // --- the last possibility: data rows inside a loop block
                if let Some(a_loop) = &mut current_loop {
                    if ls.starts_with(';') {
                        let vec = vec![read_string(ls, &mut line_iter)];
                        a_loop.add_data_row(vec);
                    }
                    else {
                        a_loop.add_data_row(split_into_strings(ls, false));
                    }
                } else {
                    panic!("Attempt to add a loop row with no loop open!");
                }
            }
        }
    }
    debug!("CIF structure loaded in: {:?}", start.elapsed());

    return data_blocks;
}

/// Returns true if the input string starts and end with a semicolon, a single or a double quote.
///
/// The input string must be trimmed!
fn is_quoted_string(input: &str) -> bool {
    let first = input.chars().nth(0).unwrap();
    if first != '\'' && first != '"' && first != ';' { return false }
    let last = input.chars().last().unwrap();
    if last != '\'' && last != '"' && last != ';' { return false }

    return true
}

fn read_string<R>(prefix: &str, line_iter: &mut Lines<R>) -> String where R: BufRead {
    let mut lines: Vec<String> = vec![prefix.to_string()];
    let mut multiline_opened = prefix.starts_with(';');
    while let Some(line) = line_iter.next() {
        match line {
            Ok(line) => {
                if line.trim().len() == 0 { continue }                  // skip empty lines
                let first_char = line.chars().nth(0).unwrap();
                if first_char != ';' && !multiline_opened {             // it's not a multiline string
                    return line;
                }
                if first_char == ';' {
                    if !multiline_opened {                              // a ';' character opens a multiline
                        lines.push(";".to_string());
                        multiline_opened = true;
                    } else {                                            // ... or closes if it has been opened
                        lines.push(";".to_string());
                        return lines.join("\n");                    // return the multiline string
                    }
                }
                // --- now the only remaining possibility is the string is a continuation of a multiline
                lines.push(line);
            }
            Err(err) => { panic!("{:?}",err); }
        }

    }

    panic!("End of CIF file reached while searching for the ';' closing a multi-line string");
}