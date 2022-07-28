use std::io::stdout;
use std::io::stderr;
use std::io::Write;
use std::path::Path;
use std::fs::{File};

/// Creates a `Writer` object.
///
/// Attempts to open a file under a given name. However, if the name is  `"stdout"` or `"stderr"`,
/// the returned `Writer` will be connected to either `stdout` or `stderr` stream, respectively.
/// Empty file name also results in writing to `stdout`.
///
/// # Arguments
///
/// * `out_fname` - file name, `"stdout"` or `"stderr"`
///
/// # Examples
///
/// ```
/// use utils::out_writer;
/// let mut to_stream = out_writer("");
/// to_stream = out_writer("stdout");
/// let mut to_file = out_writer("file.out");
/// ```
pub fn out_writer(out_fname: &str) -> Box<dyn Write>{
    match out_fname {
        "" => Box::new(stdout()) as Box<dyn Write>,
        "stdout" => Box::new(stdout()) as Box<dyn Write>,
        "stderr" => Box::new(stderr()) as Box<dyn Write>,
        _ => {
            let path = Path::new(out_fname);
            Box::new(File::options().append(true).create(true).open(&path).unwrap()) as Box<dyn Write>
        }
    }
}