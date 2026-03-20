use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use clap::{Parser, ArgAction};
use bioshell_io::{read_whitespace_delimited_values, read_delimited_values };
use bioshell_statistics::{HistogramND, into_matrix2d};

/// Command line utility for making 2D histograms
#[derive(Parser, Debug)]
struct Args {
    /// input file provides data for the histogram
    input: String,

    /// Which columns holds the X,Y data? Provide 1-based indices, otherwise the first two columns will be used
    #[arg(short = 'c', long = "columns", num_args = 2)]
    columns: Option<Vec<usize>>,

    /// Provide the observed range of data, for both dimensions
    #[arg(short = 'd', long = "data_range", num_args = 4)]
    data_range: Option<Vec<f64>>,

    /// Width for histogram bins
    #[arg(short = 'b', long = "bin_width", num_args = 2)]
    bin_width: Vec<f64>,

    // todo implement n_skip in read_delimited_values() and in read_whitespace_delimited_values()
    // /// skip the header line of the input file; note that comment lines starting with '#' are always skipped
    // #[arg(long = "skip-header", action = ArgAction::SetFalse)]
    // skip_header: bool,

    /// Output histogram file will have 3 columns: x_min y_min counts (default: hist2d.dat)
    #[arg(short = 'o', long = "output", default_value = "hist2d.dat")]
    output: String,
}

fn read_input_as_rows(path: &Path) -> io::Result<Vec<Vec<f64>>> {

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());

    let rows = match ext.as_deref() {
        Some("csv") | Some("cvs") => { read_delimited_values::<f64, _>(reader, b',')? }
        Some("tsv") => { read_delimited_values::<f64, _>(reader, b'\t')? }
        _ => read_whitespace_delimited_values::<f64, _>(reader)?
    };

    Ok(rows)
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

    // Determine 0-based column indices (default: first two columns).
    let (x_idx, y_idx) = match args.columns.as_deref() {
        Some(v) if v.len() == 2 => {
            anyhow::ensure!(v[0] >= 1 && v[1] >= 1, "Column indices must be 1-based (>= 1).");
            (v[0] - 1, v[1] - 1)
        }
        Some(_) => anyhow::bail!("--columns/-c must provide exactly two indices (x y)."),
        None => (0, 1),
    };

    let input_path = PathBuf::from(&args.input);
    let rows: Vec<Vec<f64>> = read_input_as_rows(&input_path)?;
    let bin_width = &args.bin_width;
    let mut hist = HistogramND::<2>::by_bin_widths([bin_width[0], bin_width[1]]);
    for row in &rows {
        hist.insert(&[row[x_idx], row[y_idx]]);
    }

    let (x_min, x_max, y_min, y_max);
    if let Some(v) = &args.data_range {
        (x_min, x_max, y_min, y_max) = (v[0], v[1], v[2], v[3]);
    } else {
        (x_min, x_max) = hist.data_range(0);
        (y_min, y_max) = hist.data_range(1);
    }

    let file = File::create(&args.output)?;
    let mut writer = BufWriter::new(file);
    let mat2d = into_matrix2d(&hist, x_min, x_max, y_min, y_max, 0.0);
    let [dx, dy] = hist.bin_widths();
    let mut x = x_min;
    for i in 0..mat2d.len() {
        let mut y = y_min;
        for j in 0..mat2d[i].len() {
            writeln!(writer, "{:6.2} {:6.2} {}", x, y, &mat2d[i][j])?;
            y += dy;
        }
        x += dx;
    }
    Ok(())
}
