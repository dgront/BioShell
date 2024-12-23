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
//! prefix is not considered a part of that name). That block, loaded as a [`CifData`]
//! struct, holds two key-value entries:
//! ``_name_1:value_1`` and ``_name_2:value_2``,  followed by a loop block.
//! Data items stored by a loop block are loaded into a [CifLoop] struct.
//!
//! The official specification of the CIF format can be found on
//! [this page](https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax)
//!
//! # Reading CIF files
//! A CIF file can be read using the [`read_cif_file`] function, which returns a vector of [`CifData`] structs:
//! ```
//! use bioshell_cif::{CifError, read_cif_file};
//! # fn main() -> Result<(), CifError > {
//! use bioshell_cif::CifTable;
//! let data_blocks = read_cif_file("./tests/test_data/ALA.cif")?;
//! assert_eq!(data_blocks.len(), 1);
//! // --- each dat block has a unique name:
//! assert_eq!(data_blocks[0].name(), "ALA");
//! # Ok(())
//! # }
//! ```
//!
//! Each data block contains data items as key-value pairs and may be accessed by their name.
//! You don't have to specify the full name of a data item, just make sure the provided substring is unique.
//! ```
//! # use bioshell_cif::{CifError, read_cif_file};
//! # fn main() -> Result<(), CifError > {
//! # let data_blocks = read_cif_file("./tests/test_data/ALA.cif")?;
//! let ala_block = &data_blocks[0];
//! // --- data items can be accessed by their names:
//! let id: Option<String> = ala_block.get_item("_chem_comp.id");
//! assert_eq!(id, Some("ALA".to_string()));
//! let molar_mass: Option<f64> = ala_block.get_item("_chem_comp.formula_weight");
//! # Ok(())
//! # }
//! ```
//!
//! Each loop blocks is represented as  [`CifLoop`] struct. The most convenient way however to process data stored
//! in a loop block is to use a [`CifTable`] struct.
//! ```
//! # use bioshell_cif::{CifError, read_cif_file};
//! # fn main() -> Result<(), CifError > {
//! use bioshell_cif::CifTable;
//! # let data_blocks = read_cif_file("./tests/test_data/ALA.cif")?;
//! # let ala_block = &data_blocks[0];
//! let atom_table = CifTable::new(ala_block, "_chem_comp_atom", ["atom_id","type_symbol"])?;
//! for [atom_id, atom_type] in atom_table.iter() {
//!     println!("Atom {} is {}", atom_id, atom_type);
//! }
//! # Ok(())
//! # }
//! ```
//! The example above prints out the following two columns: ``atom_id.atom_id`` and ``atom_id.type_symbol``.
//!

mod column_mapping;
mod cif_errors;
mod cif_line_iterator;

pub use column_mapping::*;
pub use cif_errors::*;
use cif_line_iterator::CifLineIterator;

use std::collections::HashMap;
use std::fmt::{Debug, Display, Formatter};
use std::io;
use std::io::{BufRead, Lines};
use std::ops::{Index, IndexMut};
use std::str::FromStr;
use std::time::Instant;
use log::{debug, info};

use bioshell_io::{open_file, split_into_strings};
use crate::CifError::{DanglingDataItem, DataValuesOutsideLoop, MisplacedDataNameInLoop, MultilineStringOutsideDataItem};

/// Returns true if a given file is in CIF format.
///
/// This function simply tests whether the first data line of a given file starts with ``data_``.
/// Otherwise, it returns ``false``. When the file can't be open returns I/O error.
pub fn is_cif_file(file_path: &str) -> io::Result<bool> {
    let reader = open_file(file_path)?;

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


/// Parses a data item into a given type or returns ItemParsingError error
///
/// # Example
/// ```
/// use bioshell_cif::CifError;
/// use bioshell_cif::CifError::ItemParsingError;
/// use bioshell_cif::parse_item_or_error;
/// fn test_macro(token: &str) -> Result<i32, CifError> {
///     let value = parse_item_or_error!(token, i32);
///     return Ok(value);
/// }
/// assert!(test_macro("1").is_ok());
/// assert!(test_macro("one").is_err());
/// ```
#[macro_export]
macro_rules! parse_item_or_error {
    ($token:expr, $type:ty) => {
        match $token.parse::<$type>() {
            Ok(val) => val,
            Err(_) => return Err( ItemParsingError {
                item: $token.to_string(), type_name: stringify!($type).to_string(), details: "".to_string(),
            }),
        }
    };
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
}

impl CifLoop {

    /// Creates an empty loop block with given columns.
    ///
    /// The newly created struct basically represents a table with named columns but with no data rows
    pub fn new(data_item_names: &[&str]) -> CifLoop {
        let cols: Vec<_> = data_item_names.iter().map(|e| e.to_string()).collect();
        return CifLoop{ column_names: cols, data_rows: vec![] };
    }

    /// Add a new column to this loop block.
    ///
    /// Adding columns is possible only before any data is inserted; once any data has been inserted,
    /// this method will panic.
    pub fn add_column(&mut self, column_name: &str) -> Result<(), CifError> {
        if self.data_rows.len() > 0 {
            return Err(MisplacedDataNameInLoop{ data_name: column_name.to_string() });
        }
        self.column_names.push(column_name.to_string());
        Ok(())
    }

    /// Add a new row of data.
    ///
    /// The provided row of data must contain the same number of entries as the number of columns
    /// in this loop-block; otherwise this method will panic.
    pub fn add_data_row(&mut self, row: Vec<String>) -> Result<(), CifError> {

        if self.last_row_complete() && row.len() == self.column_names.len() {
            self.data_rows.push(row);
            return Ok(());
        } else {
            for v in row {
                self.add_data(&v);
            }
        }

        Ok(())
    }

    /// Add a new data value to the most recent row.
    ///
    /// If the very last row contains the same number of entries as the number of columns,
    /// this method starts a new row.
    pub fn add_data(&mut self, data_value: &String) {

        if self.last_row_complete() {   // --- start a new row by inserting a new vector holding the given entry
            self.data_rows.push(vec![data_value.to_string()]);
        } else {
            let last_vec: &mut Vec<String> = self.data_rows.last_mut().unwrap();
            last_vec.push(data_value.to_string());
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

    /// Returns ``true`` if any column of this loop block contains the given substring.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_loop_example
    /// loop_
    /// _atom_site_label
    /// _atom_site_Cartn_x
    /// _atom_site_Cartn_y
    /// _atom_site_Cartn_z
    /// O1 4.154 5.699 3.026
    /// C2 5.630 5.087 4.246
    /// ";
    /// let data_blocks = read_cif_buffer(&mut BufReader::new(cif_block.as_bytes())).unwrap();
    /// # assert_eq!(data_blocks.len(), 1);
    /// # assert_eq!(data_blocks[0].name(),"loop_example");
    /// let a_loop = data_blocks[0].loop_blocks().next().unwrap();
    /// assert!(a_loop.column_name_contains("atom_site"));
    /// assert!(a_loop.column_name_contains("x"));
    /// assert!(!a_loop.column_name_contains("vector"))
    /// ```
    pub fn column_name_contains(&self, substring: &str) -> bool {
        return self.column_names.iter().any(|name| name.contains(substring));
    }

    /// Returns ``true`` if any column of this loop block starts with the given substring.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_loop_example
    /// loop_
    /// _atom_site_label
    /// _atom_site_Cartn_x
    /// _atom_site_Cartn_y
    /// _atom_site_Cartn_z
    /// O1 4.154 5.699 3.026
    /// C2 5.630 5.087 4.246
    /// ";
    /// let data_blocks = read_cif_buffer(&mut BufReader::new(cif_block.as_bytes())).unwrap();
    /// # assert_eq!(data_blocks.len(), 1);
    /// # assert_eq!(data_blocks[0].name(),"loop_example");
    /// let a_loop = data_blocks[0].loop_blocks().next().unwrap();
    /// assert_eq!(a_loop.column_name_starts_with("atom_site"), false);
    /// assert_eq!(a_loop.column_name_starts_with("x"), false);
    /// assert_eq!(a_loop.column_name_starts_with("_atom_site"), true);
    /// ```
    pub fn column_name_starts_with(&self, substring: &str) -> bool {
        return self.column_names.iter().any(|name| name.starts_with(substring));
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

    fn last_row_complete(&self) -> bool {
        if let Some(last_row) = self.data_rows.last() {
            return last_row.len() == self.column_names.len();
        }
        // if no data rows are present, we want to start one
        return true;
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
    /// let data_blocks = read_cif_buffer(&mut reader).unwrap();
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
    /// let data_blocks = read_cif_buffer(&mut reader).unwrap();
    /// assert_eq!(data_blocks.len(), 1);
    /// assert_eq!(data_blocks[0].name(),"first_block");
    /// ```
    pub fn name(&self) -> &str { &self.name }

    /// Add a given loop to this data block
    pub fn add_loop(&mut self, a_loop: CifLoop) { self.loops.push(a_loop)}

    /// Add a new data item to this data block.
    ///
    /// Note that the same results may be achieved by inserting data directly to the internal ``HashMap``
    /// returned by [``data_items_mut()``](data_items_mut()) method.
    ///
    /// # Example
    /// ```
    /// use bioshell_cif::CifData;
    /// let mut m = CifData::new("atomic_mass");
    /// m.add_item("C", "12.011");
    /// m.add_item("O", "15.999");
    /// # assert_eq!(m.data_items().len(), 2);
    /// ```
    pub fn add_item(&mut self, key: &str, data: &str) { self.data_items.insert(key.to_string(),data.to_string()); }

    /// Returns a data item value assigned to a given key.
    ///
    /// If that key can't be found in a given data block, returns ``None``.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_1AZP
    /// _entry.id                                        1AZP
    /// _refine.ls_d_res_high                            1.6
    /// ";
    /// let mut reader = BufReader::new(cif_block.as_bytes());
    /// let blocks = read_cif_buffer(&mut reader).unwrap();
    /// let data_block = &blocks[0];
    /// assert_eq!(data_block.get_item("_entry.id"), Some("1AZP".to_string()));
    /// assert_eq!(data_block.get_item("_refine.ls_d_res_high"), Some(1.6));
    /// ```
    pub fn get_item<T: FromStr>(&self, key: &str) -> Option<T> {
        self.data_items.get(key).and_then(|value| value.parse().ok())
    }

    /// Returns a data item value assigned to a given key or a default value.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// let cif_block = "data_1AZP
    /// _entry.id                                        1AZP
    /// _refine.ls_d_res_high                            1.6
    /// ";
    /// let mut reader = BufReader::new(cif_block.as_bytes());
    /// let blocks = read_cif_buffer(&mut reader).unwrap();
    /// let data_block = &blocks[0];
    /// assert_eq!(data_block.get_item_or_default("_entry.id", "1ABC".to_string()), "1AZP".to_string());
    /// assert_eq!(data_block.get_item_or_default("_refine.ls_d_res_lo", 1.8), 1.8);
    /// ```
    pub fn get_item_or_default<T: FromStr>(&self, key: &str, default_val: T) -> T where <T as FromStr>::Err: Debug {
        if let Some(value) = self.data_items.get(key) {
            value.parse().unwrap()
        } else {
            default_val
        }
    }

    /// Read access to data items of this block
    pub fn data_items(&self) -> &HashMap<String, String> { &self.data_items }

    /// Mutable access to data items of this block.
    ///
    /// This method allows to modify the set of items stored by data block
    /// # Example
    /// ```
    /// use bioshell_cif::CifData;
    /// let mut m = CifData::new("atomic_mass");
    /// m.data_items_mut().insert("O".to_string(), "15.999".to_string());
    /// # assert_eq!(m.data_items().len(), 1);
    /// ```
    pub fn data_items_mut(&mut self) -> &mut HashMap<String, String> { &mut self.data_items }

    /// Get an iterator of references to loop-blocks.
    pub fn loop_blocks(&self) -> impl DoubleEndedIterator<Item = &CifLoop> + '_ {
        self.loops.iter()
    }

    /// Finds a loop block that contains a given entry name.
    ///
    /// The provided ``data_name`` doesn't have to be a full name of an entry. Partial name is also accepted
    /// as long as it is unique.
    pub fn get_loop(&self, data_name: &str) -> Option<&CifLoop> {
        self.loops.iter().filter(|l| l.column_name_starts_with(data_name)).next()
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

impl Index<&str> for CifData {
    type Output = String;

    /// Access data items of a CifData struct with indexing operator
    ///
    /// # Example
    /// ```
    /// use bioshell_cif::CifData;
    /// let mut z = CifData::new("atomic_number");
    /// z.add_item("C", "6");
    /// assert_eq!(z["C"], "6".to_string());
    /// ```
    fn index(&self, key: &str) -> &Self::Output {
        self.data_items.get(key).expect(&format!("CifData key not found: {}", key))
    }
}

impl IndexMut<&str> for CifData {

    /// Mutable access to data items of a CifData struct with indexing operator
    ///
    ///
    /// # Example
    /// ```
    /// use bioshell_cif::CifData;
    /// let mut z = CifData::new("atomic_number");
    /// z.add_item("O", "6");
    /// z["O"] = "8".to_string();
    /// assert_eq!(z["O"], "8".to_string());
    fn index_mut(&mut self, key: &str) -> &mut Self::Output {
        self.data_items
            .entry(String::from(key))
            .or_insert(String::new())
    }
}

/// Reads a CIF-formatted file.
///
/// This function opens a file as a Returns a ``BufRead``, calls (``read_cif_buffer()``)[read_cif_buffer()]
/// and returns all the data blocks it found.
pub fn read_cif_file(input_fname: &str) -> Result<Vec<CifData>, CifError> {

    info!("Loading a CIF file: {}", input_fname);

    let reader = open_file(input_fname)?;

    return read_cif_buffer(reader);
}

enum CifLine {
    DataBlock(String),          // starts a new block, e.g. data_1AZP
    DataItem(String, String),   // e.g. _atom_site.group_PDB "A"
    LoopBlock,                  // starts a new loop block, always loop_
    MultilineString(String),    // a multi-line string
    DataName(String),           // name of a column within a loop block, e.g. _atom_site.group_PDB
    DataValues(Vec<String>),    // values of a row in a loop block, e.g. "A"
    EmptyLine,              // empty line
}

fn trim_view(s: &str) -> &str {
    let start = s.find(|c: char| !c.is_whitespace()).unwrap_or(0);
    let end = s.rfind(|c: char| !c.is_whitespace()).map_or(0, |i| i + 1);
    &s[start..end]
}

impl FromStr for CifLine {
    type Err = CifError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed_s: &str = trim_view(s);
        // --- empty string or a comment, e.g. "# this is a comment" ---
        if trimmed_s.len() == 0 || trimmed_s.starts_with("#") {
            return Ok(CifLine::EmptyLine);
        }
        // --- we start a new data block, e.g. "data_block_name" ---
        if trimmed_s.starts_with("data_") {
            return Ok(CifLine::DataBlock(trimmed_s[5..].to_string()));
        }
        // --- the name of a column in a loop block
        if trimmed_s.starts_with("_") {
            let mut data_item = split_into_strings(trimmed_s, false);
            if data_item.len() == 1 {
                return Ok(CifLine::DataName(data_item[0].to_string()));
            } else if data_item.len() == 2 {
                return Ok(CifLine::DataItem(data_item.swap_remove(0), data_item.swap_remove(0)));
            } else {
                return Err(CifError::InvalidCifLine{ line: s.to_string() });
            }
        }
        // --- a new loop block starts
        if trimmed_s.starts_with("loop_") {
            return Ok(CifLine::LoopBlock);
        }
        // --- load a multiline string
        if trimmed_s.starts_with(";") {
            return Ok(CifLine::MultilineString(trimmed_s.to_string()));
        }
        // --- load a data line split into tokens
        return Ok(CifLine::DataValues(split_into_strings(trimmed_s, false)));
    }
}


/// Reads a CIF-formatted data from a buffer.
///
/// Returns a vector of all data blocks found.
/// # Example
/// ```
/// use std::io::BufReader;
/// use bioshell_cif::read_cif_buffer;
/// let cif_block = "data_ALA
/// _chem_comp.id                                    ALA
/// _chem_comp.name                                  ALANINE
/// _chem_comp.type                                  'L-PEPTIDE LINKING'
/// _chem_comp.pdbx_type                             ATOMP
/// ";
/// let mut reader = BufReader::new(cif_block.as_bytes());
/// let data_blocks = read_cif_buffer(&mut reader);
/// assert!(data_blocks.is_ok());
/// let data_blocks = data_blocks.unwrap();
/// assert_eq!(data_blocks.len(), 1);
/// # assert_eq!(data_blocks[0].get_item("_chem_comp.id"), Some("ALA".to_string()));
/// ```
pub fn read_cif_buffer<R: BufRead>(buffer: R) -> Result<Vec<CifData>, CifError> {
    let mut data_blocks: Vec<CifData> = vec![];
    let mut current_loop: Option<CifLoop> = None;
    let mut data_item_open: Option<String> = None;

    let start = Instant::now();
    let mut line_iter = CifLineIterator::new(buffer.lines());


    while let Some(line) = line_iter.next() {
        let cif_line = CifLine::from_str(&line);
        match cif_line {

            // --- propagate errors
            Err(e) => { return Err(e); }

            // --- we have a new loop
            Ok(CifLine::LoopBlock) => {
                // --- close the previous loop if any and open a new one
                if let Some(a_loop) = current_loop {
                    add_loop_to_last_block(&mut data_blocks, a_loop)?;
                }
                current_loop = Some(CifLoop { column_names: vec![], data_rows: vec![] });
            }

            // --- we have a new data block
            Ok(CifLine::DataBlock(block_name)) => {
                // --- close the previous loop
                if let Some(a_loop) = current_loop {
                    add_loop_to_last_block(&mut data_blocks, a_loop)?;
                    current_loop = None;
                }
                data_blocks.push(CifData::new(&block_name));
            }

            // --- key-value pair
            Ok(CifLine::DataItem(key, val)) => {
                // --- close the previous loop
                if let Some(a_loop) = current_loop {
                    add_loop_to_last_block(&mut data_blocks, a_loop)?;
                    current_loop = None;
                }
                data_blocks.last_mut().unwrap().data_items_mut().insert(key, val);
            }

            Ok(CifLine::MultilineString(val)) => {
                if let Some(a_loop) = &mut current_loop {
                    a_loop.add_data(&val);
                } else {
                    if data_item_open.is_some() {
                        data_blocks.last_mut().unwrap().data_items_mut().insert(data_item_open.unwrap(), val);
                        data_item_open = None;
                    } else {
                        return Err(MultilineStringOutsideDataItem {data_value: val.to_string() });
                    }
                }
            }

            // --- data name as the only token in a line; may be a loop column or a data item
            Ok(CifLine::DataName(data_name)) => {
                if current_loop.is_some() {
                    if current_loop.as_ref().unwrap().count_rows() == 0 {
                        current_loop.as_mut().unwrap().add_column(&data_name)?;
                    } else {
                        add_loop_to_last_block(&mut data_blocks, current_loop.unwrap())?;
                        current_loop = None
                    }
                }
                if current_loop.is_none() {
                    if data_item_open.is_none() {
                        data_item_open = Some(data_name);
                    } else {
                        return Err(DanglingDataItem { data_name: data_name.to_string() });
                    }
                }
            }

            // --- a row of data values, may be a single value as well
            Ok(CifLine::DataValues(mut data_values)) => {
                if let Some(a_loop) = &mut current_loop {       // may be a loop data row
                    a_loop.add_data_row(data_values)?;
                } else {
                    if data_item_open.is_some() && data_values.len() == 1 {
                        data_blocks.last_mut().unwrap().data_items_mut().insert(data_item_open.unwrap(), data_values.swap_remove(0));
                        data_item_open = None;
                    } else {
                        return Err(DataValuesOutsideLoop { breaking_line: line.to_string() });
                    }
                }
            }
            Ok(CifLine::EmptyLine) => {}
        }
    }

    // --- close the very last loop that may be still open
    if let Some(loop_block) = current_loop {
        add_loop_to_last_block(&mut data_blocks, loop_block)?;
    }

    debug!("CIF structure loaded in: {:?}", start.elapsed());

    return Ok(data_blocks);
}

// Helper function to add a loop to the last block in the ``data_blocks`` vector.
fn add_loop_to_last_block(data_blocks: &mut Vec<CifData>, loop_block: CifLoop) -> Result<(), CifError> {
    if let Some(last_block) = data_blocks.last_mut() {
        last_block.add_loop(loop_block);
        return Ok(());
    } else { return Err(CifError::NoDataBlock); }
}

/// Returns true if a given string bears a meaningful value.
///
/// CIF format specification defines special values of '.' and '?' that represent data
/// that are inapplicable or unknown, respectively. This helper function returns ``false``
/// if a given  ``data_entry`` is one of these special characters.
///
/// # Example
/// ```
/// use bioshell_cif::entry_has_value;
/// assert!(entry_has_value("1.6"));
/// assert!(entry_has_value("A string"));
/// assert!(!entry_has_value("."));
/// assert!(!entry_has_value("?"));
/// ```
pub fn entry_has_value(data_entry: &str) -> bool { !(data_entry == "?" || data_entry == ".") }

/// Returns a value parsed from a given string whenever it's not a special symbol.
///
/// This function returns the given ``default_val`` if  [`has_value()`](entry_has_value) returns false,
/// otherwise it attempts to parse the given string into a desired type. *Note* that this function
/// doesn't check for general parsing errors, e.g. parsing ``"1.3"`` into an ``i16`` type will fail.
///
/// # Example
/// ```
/// use bioshell_cif::value_or_default;
/// assert_eq!(value_or_default::<i16>("34", 0), 34);
/// assert_eq!(value_or_default::<i16>(".", 0), 0);
/// assert_eq!(value_or_default::<i16>("?", 0), 0);
/// ```
///
pub fn value_or_default<T: FromStr>(data_entry: &str, default_val: T) -> T {
    if entry_has_value(data_entry) { data_entry.parse().ok().unwrap() } else { default_val }
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

/// Parses a string into a boolean value.
///
/// This function handles the special CIF values ``Y`` and ``N`` as well as the standard ``true`` and ``false``.
///
/// # Examples
/// ```
/// use bioshell_cif::parse_bool;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// assert_eq!(parse_bool("Y")?, true);
/// assert_eq!(parse_bool("N")?, false);
/// assert_eq!(parse_bool("true")?, true);
/// assert_eq!(parse_bool("false")?, false);
/// assert!(parse_bool("invalid").is_err());
/// # Ok(())
/// # }
/// ```
pub fn parse_bool(input: &str) -> Result<bool, CifError> {
    match input {
        "Y" | "y" => Ok(true),
        "N" | "n" => Ok(false),
        _ => input.parse::<bool>().map_err(|e| CifError::ItemParsingError {
            item: input.to_string(),
            type_name: "bool".to_string(),
            details: e.to_string(),
        }),
    }
}

