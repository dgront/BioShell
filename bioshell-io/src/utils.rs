use std::ffi::OsStr;
use std::io::{BufRead, BufReader, Error, stdout, stderr, ErrorKind, Write};
use std::path::Path;
use std::fs::{File};
use csv;
use csv::StringRecord;
use flate2::read;
use std::str::FromStr;
use std::time::Instant;
use log::debug;

/// Creates a `Writer` object.
///
/// Attempts to open a file under a given name. However, if the name is  `"stdout"` or `"stderr"`,
/// the returned `Writer` will be connected to either `stdout` or `stderr` stream, respectively.
/// Empty file name also results in writing to `stdout`.
///
/// # Arguments
/// * `out_fname` - file name, `"stdout"` or `"stderr"`
/// * `if_append` - existing file will be removed if false, otherwise the new content will be appended
///
/// # Examples
///
/// ```
/// use std::fs;
/// use bioshell_io::out_writer;
/// let mut to_stream = out_writer("", true);
/// to_stream = out_writer("stdout", true);
/// assert!(fs::metadata("stdout").is_err());
/// let mut to_file = out_writer("file.out", false);
/// assert!(fs::metadata("file.out").is_ok());
/// # fs::remove_file("file.out").expect("Can't remove a test file: file.out");
/// ```
pub fn out_writer(out_fname: &str, if_append: bool) -> Box<dyn Write>{
    match out_fname {
        "" => Box::new(stdout()) as Box<dyn Write>,
        "stdout" => Box::new(stdout()) as Box<dyn Write>,
        "stderr" => Box::new(stderr()) as Box<dyn Write>,
        _ => {
            let path = Path::new(out_fname);

            if if_append {
                let file = match File::options().append(true).write(true).create(true).open(&path) {
                    Ok(file) => file,
                    Err(e) => panic!("can't open >{:?}<, error is: {:?}", &path, e),
                } ;
                return Box::new(file) as Box<dyn Write>;
            } else {
                let file = match File::create(&path) {
                    Ok(file) => file,
                    Err(e) => panic!("can't open >{:?}<, error is: {:?}", &path, e),
                };
                return Box::new(file) as Box<dyn Write>;
            }
        }
    }
}


/// Check if all fields of the given record are not empty
fn is_record_ok(rec: &StringRecord) -> bool {

    let mut flag = true;
    for e in rec {
        if e.len() == 0 {
            flag = false;
            break;
        }
    }
    return flag;
}

/// Reads values from a file that are delimited with a given character.
///
/// The function can handle both tab-separated and comma-separated files.
///
/// # Examples
///
/// Read a tab-separated file:
/// ```
/// use std::io::BufReader;
/// # use std::io;
/// use bioshell_io::{open_file, read_delimited_values};
/// # fn main() -> Result<(), io::Error> {
/// let txt_f64 = "1.0\t2.0\t3.0\t4.0
/// 5.0\t6.0\t7.0\t8.0
/// 9.0\t10.0\t11.0\t12.0
/// ";
/// let data_f64: Vec<Vec<f64>> = read_delimited_values(BufReader::new(txt_f64.as_bytes()), b'\t').unwrap();
/// assert_eq!(data_f64.len(), 3);
/// assert_eq!(data_f64[0].len(), 4);
/// let buffer = open_file("tests/test_files/string.tsv")?;
/// let data_str: Vec<Vec<String>> = read_delimited_values(buffer, b'\t')?;
/// # assert_eq!(data_str.len(), 3);
/// # assert_eq!(data_str[0].len(), 2);
/// # Ok(())
/// # }
/// ```
///
/// Read a comma-separated file:
/// ```
/// # use std::io;
/// # fn main() -> Result<(), io::Error> {
/// use bioshell_io::{open_file, read_delimited_values};
/// let reader = open_file("tests/test_files/f64.csv")?;
/// let data_f64: Vec<Vec<f64>> = read_delimited_values(reader, b',')?;
/// # assert_eq!(data_f64.len(), 2);
/// # assert_eq!(data_f64[1].len(), 3);
/// # Ok(())
/// }
/// ```
pub fn read_delimited_values<T: FromStr+Clone, R: BufRead>(reader: R, delimiter:u8) -> Result<Vec<Vec<T>>, Error> {

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .comment(Some(b'#'))
        .from_reader(reader);

    let mut table: TableOfRows<T> = TableOfRows{ table: Table::new_empty() };
    for record in rdr.records() {
        if let Ok(r) = &record {
            if !is_record_ok(r) { continue; }
            table.insert_row(r.into_iter())?;
        }
    }

    return Ok(table.table.data);
}

/// Reads a file that is delimited with a given character and returns data loaded column-wise
///
/// The function can handle both tab-separated and comma-separated files.
///
/// # Examples
///
/// Read a tab-separated file:
/// ```
/// use std::io::BufReader;
/// # use std::io;
/// use bioshell_io::{open_file, read_delimited_columns};
/// # fn main() -> Result<(), io::Error> {
/// let txt_f64 = "1.0\t2.0\t3.0\t4.0
/// 5.0\t6.0\t7.0\t8.0
/// 9.0\t10.0\t11.0\t12.0
/// ";
/// let data_f64: Vec<Vec<f64>> = read_delimited_columns(BufReader::new(txt_f64.as_bytes()), b'\t').unwrap();
/// assert_eq!(data_f64.len(), 4);
/// assert_eq!(data_f64[0].len(), 3);
/// let buffer = open_file("tests/test_files/string.tsv")?;
/// let data_str: Vec<Vec<String>> = read_delimited_columns(buffer, b'\t')?;
/// # assert_eq!(data_str.len(), 2);
/// # assert_eq!(data_str[0].len(), 3);
/// # Ok(())
/// # }
/// ```
pub fn read_delimited_columns<T: FromStr+Clone, R: BufRead>(reader: R, delimiter:u8) -> Result<Vec<Vec<T>>, Error> {

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .comment(Some(b'#'))
        .from_reader(reader);

    let mut table: TableOfColumns<T> = TableOfColumns { table: Table::new_empty() };
    for record in rdr.records() {
        if let Ok(r) = &record {
            if !is_record_ok(r) { continue; }
            table.insert_row(r.into_iter())?;
        }
    }

    return Ok(table.table.data);
}


/// Reads whitespace-separated data from a `BufReader` and stores it in a `Vec<Vec<T>>`.
///
/// The function assumes that the number of columns is constant across all rows.
pub fn read_whitespace_delimited_values<T, R>(reader: R) -> Result<Vec<Vec<T>>, Error>
where T: FromStr+Clone,R: BufRead,
{
    let mut table: TableOfRows<T> = TableOfRows{ table: Table::new_empty() };
    for line in reader.lines() {
        let line = line?;  // Read the line (handle I/O errors)
        if line.is_empty() || line.starts_with('#') { // Skip empty lines and comments
            continue;
        }
        let fields_iter = line.split_whitespace();
        table.insert_row(fields_iter)?;
    }

    Ok(table.table.data)
}

/// Reads whitespace-separated data from a `BufReader` and stores it column-wise  in a `Vec<Vec<T>>`.
///
/// The function assumes that the number of columns is constant across all rows.
pub fn read_whitespace_delimited_columns<T, R>(reader: R) -> Result<Vec<Vec<T>>, Error>
where T: FromStr+Clone,R: BufRead,
{
    let mut table: TableOfColumns<T> = TableOfColumns { table: Table::new_empty() };
    let start = Instant::now();

    for line in reader.lines() {
        let line = line?;  // Read the line (handle I/O errors)
        if line.is_empty() || line.starts_with('#') { // Skip empty lines and comments
            continue;
        }
        table.insert_row(line.split_whitespace())?;
    }
    debug!("text file loaded in: {:?}", start.elapsed());

    Ok(table.table.data)
}

struct Table<T:FromStr> { data: Vec<Vec<T>>, }

trait InsertRow<T: FromStr> {
    fn insert_row<'a>(&mut self, row: impl Iterator<Item=&'a str>) -> Result<(), Error>;
}

impl<T: FromStr+Clone> Table<T> {
    pub fn new_empty() -> Self { Table { data: vec![] } }
    pub fn new_allocated(n_rows: usize, n_columns: usize, val: T) -> Self {
        Table { data:  vec![vec![val; n_rows]; n_columns]}
    }
}

struct TableOfRows<T:FromStr> { table: Table<T>, }

impl<T: FromStr> InsertRow<T> for TableOfRows<T> {
    fn insert_row<'a>(&mut self, row: impl Iterator<Item=&'a str>) -> Result<(), Error> {
        let row: Vec<T> = row.into_iter()
            .map(|field| {
                field.parse::<T>().map_err(|_e| {
                    Error::new(ErrorKind::InvalidData, format!("Failed to parse field: {}", field))
                })
            }).collect::<Result<Vec<T>, Error>>()?;

        self.table.data.push(row);
        Ok(())
    }
}

struct TableOfColumns<T:FromStr+Clone> { table: Table<T>, }


impl<T: FromStr+Clone> InsertRow<T> for TableOfColumns<T> {
    fn insert_row<'a>(&mut self, row: impl Iterator<Item=&'a str>) -> Result<(), Error> {

        let row: Vec<T> = row.into_iter()
            .map(|field| {
                field.parse::<T>().map_err(|_e| {
                    Error::new(ErrorKind::InvalidData, format!("Failed to parse field: {}", field))
                })
            }).collect::<Result<Vec<T>, Error>>()?;

        if self.table.data.len() == 0 {
            for v in row.into_iter() { self.table.data.push(vec![v]); }
        } else { for (i, v) in row.into_iter().enumerate() { self.table.data[i].push(v); } }

        Ok(())
    }
}

/// Opens a file for reading.
///
/// This function can open a regular file or a gzipped one, as determined by the extension
/// of the input file name. A boxed reader to the content is returned.
///
/// # Examples
/// ```
/// use bioshell_io::open_file;
/// # use std::io;
/// # fn main() -> Result<(), io::Error> {
/// let reader = open_file("tests/test_files/f64.csv")?;
/// let reader_gzipped = open_file("tests/test_files/f64.csv.gz")?;
/// # Ok(())
/// # }
/// ```
pub fn open_file(filename: &str) -> Result<Box<dyn BufRead>, Error> {
    if filename.len() == 0 {
        panic!("\nCouldn't open file - file name is an empty string!");
    }
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => {return Err(why)},
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Ok(Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(128 * 1024, file)))
    }
}

/// Counts the number of rows in a text file.
pub fn count_rows(file_path: &str) -> Result<usize, Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let row_count = reader.lines().count();

    Ok(row_count)
}