use std::ffi::OsStr;
use std::io::{BufRead, BufReader, stdout};
use std::io::stderr;
use std::io::Write;
use std::path::Path;
use std::fs::{File};
use csv;
use csv::StringRecord;
use flate2::read;

/// Check whether a `Writer` object created based on a given string would actually write to screen.
///
/// Returns true if the name is  `"stdout"` `"stderr"` or empty.
///
/// # Arguments
///
/// * `out_fname` - file name, `"stdout"` or `"stderr"`
///
/// # Examples
///
/// ```
/// use bioshell_io::writes_to_screen;
/// assert!(writes_to_screen(""));
/// assert!(writes_to_screen("stderr"));
/// assert!(writes_to_screen("stdout"));
/// assert!(!writes_to_screen("file.txt"));
/// ```
pub fn writes_to_screen(out_fname: &str) -> bool {
    match out_fname {
        "" =>true,
        "stdout" => true,
        "stderr" => true,
        _ => false
    }
}

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
/// use bioshell_io::out_writer;
/// let mut to_stream = out_writer("", true);
/// to_stream = out_writer("stdout", true);
/// let mut to_file = out_writer("file.out", false);
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

/// Reads real values from a file in the tab-separated format
pub fn read_tsv(fname: &str) -> Vec<Vec<f64>> { read_csv_tsv(fname, b'\t') }

/// Reads real values from a file in the coma-separated format
pub fn read_csv(fname: &str) -> Vec<Vec<f64>> { read_csv_tsv(fname, b',') }

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

fn read_csv_tsv(fname:&str, delimiter:u8) -> Vec<Vec<f64>> {

    // --- this BioShell utility handles all I/O expceptions
    let reader = open_file(fname);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .from_reader(reader);
    let mut data : Vec<Vec<f64>> = Vec::new();
    for record in rdr.records() {
        if let Ok(r) = &record {
            if !is_record_ok(r) { continue; }
            let row : Vec<f64> = r.iter().map(|e| {
                match e.parse::<f64>(){
                    Ok(v) => v,
                    Err(_err) => panic!("Problem while parsing a float value: {}\nThe last record was: {:?}", e, &record),
                }
            }).collect();
            data.push(row);
        }
    }

    return data;
}

/// Opens a file for reading.
///
/// This function can open a regular file or a gzipped one, as determined by the extension
/// of the input file name. A boxed reader to the content is returned.
pub fn open_file(filename: &str) -> Box<dyn BufRead> {
    if filename.len() == 0 {
        panic!("\nCouldn't open file - file name is an empty string!");
    }
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("\nCouldn't open file '{}': {}", path.display(), why),
        Ok(file) => file,
    };

    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(
            128 * 1024,
            read::GzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}

/// Splits a string by a whitespace into strings that can contain a whitespace within
pub fn split_into_strings(s: &str) -> Vec<String> {
    let mut wrapped = false;

    s.split(|c| {
        if c == '"' {
            wrapped = !wrapped;
        }
        c == ' ' && !wrapped
    })
        // .map(|s| s.replace("\"", "")) // remove the quotation marks from the sub-strings
        .map(|s|s.to_string())
        .filter(|s| s.len() > 0)
        .collect()
}