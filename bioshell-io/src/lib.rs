//! Utility functions to facilitate I/O operations for bioshell crates
//!
//! A few functions that have been frequently used by bioshell crates were refactored
//! and gathered within ``bioshell-io`` crate. While the set of these functions will
//! be most likely growing, currently the crate allows for:
//!
//! # Opening an input stream, which might be gzip'ed
//!
//! The [open_file()] function opens a file of a given name. If that file name
//! has ``.gz`` suffix, the returned ``BufRead`` is automatically uncompressed:
//!
//! ```
//! use bioshell_io::open_file;
//! # use std::io;
//! # fn main() -> Result<(), io::Error> {
//! let reader = open_file("tests/test_files/f64.csv")?;
//! let reader_gzipped = open_file("tests/test_files/f64.csv.gz")?;
//! # Ok(())
//! # }
//! ```
//!
//! # Unified opening an output stream.
//!
//! [out_writer()] opens a file for writing. If the given file name is ``"stdout"`` or ``"stderr"``, writes to the appropriate
//! stream rather than to a file:
//!
//! ```
//! use std::fs;
//! use bioshell_io::out_writer;
//! // This will print on stdout
//! let mut to_stream = out_writer("", true);
//! // This will also print on stdout
//! to_stream = out_writer("stdout", true);
//! // "stdout" file should not exist
//! assert!(fs::metadata("stdout").is_err());
//! // now let's open a regular file for writing
//! let mut to_file = out_writer("file.out", false);
//! assert!(fs::metadata("file.out").is_ok());
//! # fs::remove_file("file.out").expect("Can't remove a test file: file.out");
//! ```
//!
//! # Reading ``.csv`` and ``.tsv`` files
//!
//! Actually, bioshell utilizes ``csv`` crate to read ``.csv`` and ``.tsv`` files. The extra job
//! that [read_delimited_values()] function do is automated parsing to a table of statically defined type, e.g. ``f64``.
//! The [read_delimited_columns()] function stores values column-wise.
//!
//! ```
//! # use std::io;
//! # fn main() -> Result<(), io::Error> {
//! use bioshell_io::{open_file, read_delimited_columns, read_delimited_values};
//! let reader = open_file("tests/test_files/f64.csv")?;
//! let data_f64: Vec<Vec<f64>> = read_delimited_values(reader, b',')?;
//! # assert_eq!(data_f64.len(), 2);
//! # assert_eq!(data_f64[1].len(), 3);
//! let reader = open_file("tests/test_files/f64.csv")?;
//! let columns_f64: Vec<Vec<f64>> = read_delimited_columns(reader, b',')?;
//! # assert_eq!(columns_f64.len(), 3);
//! # assert_eq!(columns_f64[1].len(), 2);
//! # Ok(())
//! # }
//! ```
//!
//! # Reading white-space delimited files
//!
//! Flat text tables of values delimited by white-space characters may be loaded
//! both row-wise with  [read_whitespace_delimited_values()] and column-wise with [read_whitespace_delimited_columns()]
//!
//! ```
//! # use std::io;
//! # fn main() -> Result<(), io::Error> {
//! use bioshell_io::{open_file, read_whitespace_delimited_columns, read_whitespace_delimited_values};
//! let reader = open_file("tests/test_files/numbers.txt")?;
//! let data_f64: Vec<Vec<f64>> = read_whitespace_delimited_values(reader)?;
//! # assert_eq!(data_f64.len(), 2);
//! # assert_eq!(data_f64[1].len(), 3);
//! let reader = open_file("tests/test_files/numbers.txt")?;
//! let columns_f32: Vec<Vec<f32>> = read_whitespace_delimited_columns(reader)?;
//! # assert_eq!(columns_f32.len(), 3);
//! # assert_eq!(columns_f32[1].len(), 2);
//! # Ok(())
//! # }
//! ```

//!
//! # Splitting a string into tokens by whitespaces
//!
//! Sure, rust provides such a functionality. The [`str::split_whitespace()`](str::split_whitespace) method however can't
//! split in complex scenarios, such as: ``"first second 'third token'"``. Even more complex example
//! is given below:
//! ```
//! use bioshell_io::split_into_strings;
//! let tokens = split_into_strings("A   'RNA linking'       y \"ADENOSINE-5'-MONOPHOSPHATE\" ? 'C10 H14 N5 O7 P' 347.221", false);
//! assert_eq!(tokens.len(), 7);
//! assert_eq!(tokens[3], "\"ADENOSINE-5'-MONOPHOSPHATE\"".to_string());
//!
//! ```
//!
#![allow(clippy::needless_return)]
mod utils;
mod split_into_strings;

pub use utils::*;
pub use split_into_strings::*;
