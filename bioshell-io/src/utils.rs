use std::{env, fs, io};
use std::ffi::OsStr;
use std::io::{BufRead, BufReader, Error, stdout, stderr, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::fs::{File};
use csv;
use csv::StringRecord;
use flate2::read;
use std::str::FromStr;
use std::time::Instant;
use log::{debug, info};

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

// -------------------------------------------------------------------------------------------------------------------
// Private structs used by read_whitespace_delimited_columns(), read_whitespace_delimited_values() and read_delimited_columns()
// -------------------------------------------------------------------------------------------------------------------

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

/// Attempts to create a file at the specified path to check if it can be written.
///
/// This function tries to create (or truncate) the file at the given path.
/// If the file is successfully created, it is immediately deleted to avoid side effects.
///
/// # Arguments
///
/// * `path` - A reference to a path that specifies where to attempt file creation.
///
/// # Returns
///
/// * `true` if the file can be created for writing.
/// * `false` if the file cannot be created (e.g., due to permission issues, invalid path, etc.).
///
/// # Notes
///
/// * If the file already exists, it **will be truncated** during this check.
/// * This operation requires appropriate write permissions for the target directory.
pub fn can_create_file<P: AsRef<Path>>(path: P) -> bool {
    match File::create(&path) {
        Ok(_) => {
            let _ = std::fs::remove_file(&path);    // --- File successfully created, now remove it
            true
        }
        Err(_) => false,
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
pub fn open_file<P: AsRef<Path>>(file_path: P) -> Result<Box<dyn BufRead>, Error> {

    if file_path.as_ref().as_os_str().is_empty() {
        panic!("\nCouldn't open file - file name is an empty string!");
    }
    let file = match File::open(&file_path) {
        Err(why) => {return Err(why)},
        Ok(file) => file,
    };

    if file_path.as_ref().extension() == Some(OsStr::new("gz")) {
        Ok(Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::with_capacity(128 * 1024, file)))
    }
}

/// Returns the full extension of a file path.
///
/// Unlike the `path.extension()` function, this function returns the full extension, including all its components.
///
/// ```
/// use std::path::PathBuf;
/// use bioshell_io::full_extension;
/// assert_eq!(&full_extension(&PathBuf::from("test.txt")).unwrap(), "txt");
/// assert_eq!(&full_extension(&PathBuf::from("file.tar.gz")).unwrap(), "tar.gz");
/// ```
pub fn full_extension(path: &PathBuf) -> Option<String> {
    let file_name = path.file_name()?.to_str()?;
    let parts: Vec<&str> = file_name.split('.').collect();

    if parts.len() > 1 {
        Some(parts[1..].join(".")) // Join all parts after the first dot
    } else {
        None
    }
}

/// Counts the number of rows in a text file.
pub fn count_rows(file_path: &str) -> Result<usize, Error> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let row_count = reader.lines().count();

    Ok(row_count)
}

/// Try to find the main BioShell v.4 folder.
///
/// This function first checks if the environment variable `BIOSHELL4_PATH` is set and points to a valid folder.
/// If not, it checks a few most obvious locations (e.g. the current folder).
///
/// The BioShell v.4 path is necessary to find folders with data used internally by BioShell crates,
/// e.g. the `bioshell-pdb` crate uses substitution matrices to score a sequence alignment
/// or `.cif` files with monomer definitions
///
/// ```
/// let bioshell_main_dir = bioshell_io::find_bioshell_path();
/// assert!(bioshell_main_dir.is_some());
/// let bioshell_main_dir = bioshell_main_dir.unwrap();
/// ```
pub fn find_bioshell_path() -> Option<PathBuf> {
    // Check if the environment variable is set and points to a valid folder.
    if let Ok(env_path) = env::var("BIOSHELL4_PATH") {
        let env_path = PathBuf::from(env_path);
        if env_path.exists() && env_path.is_dir() {
            info!("Main BioShell v.4 path located as {}", env_path.to_str().unwrap());
            return Some(env_path);
        }
    }

    // Check if the current path is the main bioshell folder
    if let Ok(current_dir) = env::current_dir() {
        let data_path = current_dir.join("README.md");
        if data_path.exists() {
            info!("Main BioShell v.4 path located as {}", current_dir.to_str().unwrap());
            return Some(current_dir);
        }

        // Check one folder up
        if let Some(up_dir) = current_dir.parent() {
            let data_path = up_dir.join("README.md");
            if data_path.exists() {
                info!("Main BioShell v.4 path located as {}", up_dir.to_str().unwrap());
                return Some(up_dir.to_path_buf());
            }
        }
    }

    None
}

/// Recursively finds files in a directory matching an extension specified by a regex.
///
/// # Arguments
/// - `dir`: The directory to search.
/// - `extension`: The file extension regex to match (e.g., "rs" for Rust files).
///
/// # Returns
/// A `Result<Vec<PathBuf>, io::Error>` containing the paths of matching files or an error.
///
/// # Example
/// ```
/// use std::path::Path;
/// use bioshell_io::glob;
/// let toml_files = glob(Path::new("./"), "toml").unwrap();
/// assert_eq!(toml_files.len(), 1);
/// let csv_gz_files = glob(Path::new("./"), r"csv\.gz").unwrap();
/// assert_eq!(csv_gz_files.len(), 1);
/// let csv_gz_files = glob(Path::new("./"), r"toml$|csv\.gz$").unwrap();
/// assert_eq!(csv_gz_files.len(), 2);
/// ```
pub fn glob(dir: &Path, extension_regex: &str) -> Result<Vec<PathBuf>, Error> {
    let mut results = Vec::new();

    // Compile the regular expression
    let regex = regex::Regex::from_str(extension_regex).map_err(|_| {
        Error::new(io::ErrorKind::InvalidInput, "Failed to parse the regex pattern")
    })?;

    if dir.is_dir() {
        // Read the directory entries
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();

            // If the path is a directory, recurse into it
            if path.is_dir() {
                results.extend(glob(&path, extension_regex)?);
            } else if let Some(file_name) = path.file_name() {
                // Check if the file name matches the regex
                if regex.is_match(file_name.to_string_lossy().as_ref()) {
                    results.push(path);
                }
            }
        }
    }

    Ok(results)
}

/// Converts a Markdown string to plain text.
///
/// The function uses [`pulldown-cmark`](https://docs.rs/pulldown_cmark) crate to parse a markdown file.
/// It attempts to preserve the formatting of the Markdown text. Handles newlines, lists, and other formatting elements.
///
/// # Example
/// ```
/// let markdown = "# Heading
/// - List item 1
/// - List item 2
/// ";
/// let text = bioshell_io::markdown_to_text(markdown);
/// assert_eq!(text, "Heading\n- List item 1\n- List item 2\n");
/// ```
pub fn markdown_to_text(md: &str) -> String {

    fn new_list_item(lists_started: &mut Vec<Option<u64>>) -> Option<u64> {
        if lists_started.is_empty() { return None; }
        if let Some(last) = lists_started.last_mut() {
            if let Some(n) = last {
                let current = *n;
                *last = Some(current + 1);
                return Some(current);
            }
        }

        None
    }

    let parser = pulldown_cmark::Parser::new(md);
    let mut plain_text = String::new();
    let mut lists_started: Vec<Option<u64>> = vec![]; // Tracks started lists
    for event in parser {
        match event {
            pulldown_cmark::Event::Start(pulldown_cmark::Tag::CodeBlock(_)) => {
                plain_text.push('\n'); // Start a new line before the code block
            }
            pulldown_cmark::Event::End(pulldown_cmark::Tag::CodeBlock(_)) => {
                plain_text.push('\n'); // End the code block with a newline
            }
            pulldown_cmark::Event::Text(text) => {
                plain_text.push_str(&text);
            }
            pulldown_cmark::Event::Start(pulldown_cmark::Tag::List(Some(start))) => {
                lists_started.push(Some(start)); // Start ordered list at given number
            }
            pulldown_cmark::Event::Start(pulldown_cmark::Tag::List(None)) => {
                lists_started.push(None); // Start unordered list
            }
            pulldown_cmark::Event::Start(pulldown_cmark::Tag::Item) => {

                if let Some(ref mut num) = new_list_item(&mut lists_started) {
                    plain_text.push_str(&format!("\n{}. ", num)); // Add numbered list item
                } else {
                    plain_text.push_str("\n- "); // Bullet point for unordered lists
                }
            }
            pulldown_cmark::Event::End(pulldown_cmark::Tag::List(_)) => {
                lists_started.pop(); // End list
            }
            pulldown_cmark::Event::SoftBreak => {
                plain_text.push('\n'); // Soft line breaks
            }
            pulldown_cmark::Event::HardBreak => {
                plain_text.push('\n'); // Explicit line breaks
            }
            pulldown_cmark::Event::End(pulldown_cmark::Tag::Paragraph) => {
                plain_text.push('\n'); // Preserve paragraph separation
            }
            _ => {}
        }
    }
    if !plain_text.ends_with('\n') { plain_text.push('\n'); }

    plain_text.to_string()
}

/// Sanitizes a string for safe use as a file name on Windows, Linux, and macOS.
///
/// Replaces or removes characters that are not allowed in file names on most platforms.
/// This includes characters like `<>:"/\\|?*` and control characters.
///
/// # Example
/// ```
/// use bioshell_io::sanitize_filename;
/// let name = "RefSeq|XP_123456.1/SwissProt:Q9NQX5.fasta";
/// assert_eq!(sanitize_filename(name), "RefSeq_XP_123456.1_SwissProt_Q9NQX5.fasta");
/// ```
pub fn sanitize_filename<S: AsRef<str>>(name: S) -> String {
    name.as_ref()
        .chars()
        .filter_map(|c| match c {
            // Disallowed on Windows and/or POSIX systems
            '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => Some('_'),
            // Control characters (ASCII 0–31) are unsafe
            c if c.is_control() => None,
            // Otherwise safe
            c => Some(c),
        })
        .collect()
}
