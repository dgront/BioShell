use std::ffi::OsStr;
use std::io::{BufRead, BufReader, Error, stdout, ErrorKind};
use std::io::stderr;
use std::io::Write;
use std::path::Path;
use std::fs::{File};
use csv;
use csv::StringRecord;
use flate2::read;


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

/// Reads values from a file in the tab-separated format
///
/// # Example
/// ```
/// use std::io::BufReader;
/// # use std::io;
/// use bioshell_io::{open_file, read_tsv};
/// # fn main() -> Result<(), io::Error> {
/// let txt_f64 = "1.0\t2.0\t3.0\t4.0
/// 5.0\t6.0\t7.0\t8.0
/// 9.0\t10.0\t11.0\t12.0
/// ";
/// let data_f64: Vec<Vec<f64>> = read_tsv(BufReader::new(txt_f64.as_bytes())).unwrap();
/// assert_eq!(data_f64.len(), 3);
/// assert_eq!(data_f64[0].len(), 4);
/// let buffer = open_file("tests/test_files/string.tsv")?;
/// let data_str: Vec<Vec<String>> = read_tsv(buffer)?;
/// # assert_eq!(data_str.len(), 2);
/// # assert_eq!(data_str[0].len(), 2);
/// # Ok(())
/// # }
/// ```
pub fn read_tsv<T: std::str::FromStr, R: BufRead>(reader: R) -> Result<Vec<Vec<T>>, Error> { read_csv_tsv(reader, b'\t') }

/// Reads values from a file in the coma-separated format
///
/// This function works as [read_tsv()], just with another delimiter
///
/// # Example
/// ```
/// # use std::io;
/// # fn main() -> Result<(), io::Error> {
/// use bioshell_io::{open_file, read_csv};
/// let reader = open_file("tests/test_files/f64.csv")?;
/// let data_f64: Vec<Vec<f64>> = read_csv(reader)?;
/// # assert_eq!(data_f64.len(), 2);
/// # assert_eq!(data_f64[1].len(), 3);
/// # Ok(())
/// }
/// ```
pub fn read_csv<T: std::str::FromStr, R: BufRead>(reader: R) -> Result<Vec<Vec<T>>, Error> { read_csv_tsv(reader, b',') }

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

fn read_csv_tsv<T: std::str::FromStr, R: BufRead>(reader: R, delimiter:u8) -> Result<Vec<Vec<T>>, Error> {

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(delimiter)
        .from_reader(reader);

    let mut data : Vec<Vec<T>> = Vec::new();
    for record in rdr.records() {
        if let Ok(r) = &record {
            if !is_record_ok(r) { continue; }

            let row: Result<Vec<T>, _> = r.iter().map(|e| {
                e.parse::<T>()
            }).collect();

            let row = match row {
                Ok(values) => values,
                Err(_err) => {
                    return Err(Error::new(ErrorKind::Other, format!("Problem while parsing a float value; the last record was: {:?}", &record)));
                }
            };

            data.push(row);
        }
    }

    return Ok(data);
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


/// Splits a string by a whitespace into strings that can contain a whitespace within
///
/// Unlike the [`str::split_whitespace()`](str::split_whitespace) from the Rust standard library, this function
/// takes quotation marks (both single and double) into account; a quoted string, that may contain white space
/// characters, is returned as a single token. The function also allows for nested quotes.
/// The second parameter of the function determines whether the quotation mark are removed
/// from a substring token or not.
///
/// # Examples
///
/// ```rust
/// use bioshell_io::split_into_strings;
/// let tokens_easy = split_into_strings("The quick brown fox jumps over the lazy dog", false);
/// assert_eq!(tokens_easy.len(), 9);
/// let tokens_quoted = split_into_strings("The 'quick brown fox' jumps over the 'lazy dog'", false);
/// assert_eq!(tokens_quoted.len(), 6);
/// assert_eq!(tokens_quoted[1], "'quick brown fox'");
/// let tokens_quoted = split_into_strings("The 'quick brown fox' jumps over the 'lazy dog'", true);
/// assert_eq!(tokens_quoted[1], "quick brown fox");
/// let tokens_tricky = split_into_strings("O \"O5'\" \"O5'\"", false);
/// assert_eq!(tokens_tricky.len(), 3);
/// let tokens_tricky = split_into_strings("O \"O5'\" \"O1\"", false);
/// assert_eq!(tokens_tricky.len(), 3);
/// let cif_tokens = split_into_strings("A   'RNA linking'       y \"ADENOSINE-5'-MONOPHOSPHATE\" ? 'C10 H14 N5 O7 P' 347.221", false);
/// assert_eq!(cif_tokens.len(), 7);
/// assert_eq!(cif_tokens[3], "\"ADENOSINE-5'-MONOPHOSPHATE\"".to_string());
/// ```
pub fn split_into_strings(s: &str, if_remove_quotes: bool) -> Vec<String> {

    let mut tokens: Vec<String> = Vec::new();
    let mut iter = s.split_whitespace();
    tokens.push(iter.next().unwrap().to_string());
    for word in iter {
        let mut token = tokens.last_mut().unwrap();
        match quote_style(token) {
            QuoteStyle::None | QuoteStyle::End(_) | QuoteStyle::Both(_) => { tokens.push(word.to_string()); }
            QuoteStyle::Begin(c) => {
                token.push_str(" ");
                token.push_str(word);
            }
        }
    }

    return tokens;
}

#[derive(PartialEq, Debug)]
enum QuoteStyle {
    None,
    Begin(char),
    End(char),
    Both(char),
}
fn quote_style(token: &str) -> QuoteStyle {
    let first_char = token.chars().next();
    let last_char = token.chars().rev().next();

    match (first_char, last_char) {
        (Some(f), Some(l)) if f == l && is_quote_char(f) => QuoteStyle::Both(f),
        (Some(f), _) if is_quote_char(f) => QuoteStyle::Begin(f),
        (_, Some(l)) if is_quote_char(l) => QuoteStyle::End(l),
        _ => QuoteStyle::None,
    }
}

/// Checks if a character is a common quotation mark
fn is_quote_char(c: char) -> bool { c == '"' || c == '\'' || c == '“' || c == '”' || c == '‘' || c == '’' }

#[cfg(test)]
mod tests {
    use crate::utils::{quote_style, QuoteStyle};

    #[test]
    fn test_quote_style() {
        let examples = [
            ("\"Hello\"", QuoteStyle::Both('"')),
            ("'World", QuoteStyle::Begin('\'')),
            ("Rust'", QuoteStyle::End('\'')),
            ("No quotes", QuoteStyle::None),
            ("\"Fancy\"", QuoteStyle::Both('"')),
            ("'Smart'", QuoteStyle::Both('\'')),
            ("'Different\"", QuoteStyle::Begin('\'')),
        ];

        for (input, expected) in examples {
            assert_eq!(quote_style(input), expected, "Failed on input: {}", input);
        }
    }
}
