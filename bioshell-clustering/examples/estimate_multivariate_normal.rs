use std::io::Error;
use clap::{Parser};

use bioshell_statistics::{Estimable, MultiNormalDistribution};
use bioshell_io::{open_file, read_delimited_values};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Estimates multivariate normal distribution parameters: an average vector and a covariance matrix
/// from a file with observed sample
/// say estimate_multivariate_normal -h to see options
struct Args {
    /// file with N-dimensional input observations: N columns of real values
    #[clap(short, long, short='i', required = true)]
    infile: String,
}

fn main() -> Result<(), Error> {

    let args = Args::parse();

    let fname = &args.infile;
    let reader = open_file(fname).expect(&format!("Can't open {} file!", fname));
    let sample = read_delimited_values(reader, b'\t')?;
    let n_dim: usize = sample[0].len();
    let n_data = sample.len();
    println!("{} rows loaded, data dimension is {}", n_data, n_dim);

    // ---------- Distribution to be inferred
    let mut normal = MultiNormalDistribution::new(n_dim);
    normal.estimate(&sample);
    println!("{}", &normal);

    Ok(())
}