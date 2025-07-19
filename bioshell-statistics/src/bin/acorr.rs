use std::env;
use std::time::Instant;
use clap::Parser;
use log::info;
use bioshell_io::{open_file, read_delimited_columns, read_whitespace_delimited_columns};
use bioshell_statistics::autocorrelate_vectors;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Command line tool to calculate autocorrelation of a time series
/// say pdb_tool -h to see options
struct Args {
    /// input data in TSV format
    #[clap(long, value_name = "FILE")]
    tsv: Option<String>,
    /// input data in whitespace separated columns
    #[clap(short='i', long, value_name = "FILE")]
    infile: Option<String>,
    /// the input is 3D vector data
    #[clap(long)]
    vector: bool,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() {
    // Parse the command-line arguments
    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    if args.verbose {
        env::set_var("RUST_LOG", "debug");
    }
    env_logger::init();

    let mut columns: Vec<Vec<f64>> = Vec::new();
    if let Some(infile) = &args.infile {
        let reader = open_file(infile).expect(&format!("Can't open {} file!", infile));
        columns = read_whitespace_delimited_columns(reader).expect("Can't parse a flat text file!");
    }
    if let Some(tsv) = &args.tsv {
        let reader = open_file(tsv).expect(&format!("Can't open {} file!", tsv));
        columns = read_delimited_columns(reader, b'\t').expect("Can't parse a .tsv file!");
    }
    if columns.len() < 1 {
        panic!("No data provided!");
    }


    // Check if --vector flag is provided
    if args.vector {
        info!("Processing as 3D vector data...");
        let start = Instant::now();
        let mut results = Vec::new();
        let num_triplets = columns.len() / 3;
        for i in 0..num_triplets {
            results.push(autocorrelate_vectors(&columns[3 * i], &columns[3 * i + 1], &columns[3 * i + 2]));
        }
        info!("{} autocorrelation functions of 3D vectors computed in: {:?}", num_triplets, start.elapsed());
        for lag in 0..(results[0].len() as f64 / 1000.0) as usize {
            print!("{:<6}", lag); // Print lag time in the first column
            for result in &results {
                print!("{:<9.6}", result[lag]); // Print each autocorrelation result for the given lag
            }
            println!(); // New line after each lag time row
        }
    } else {
        info!("Processing as scalar data...");
    }

}
