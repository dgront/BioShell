use std::env;
use std::iter::zip;
use clap::{Parser};
use log::{info};

use bioshell_clustering::kmeans::KMeans;
use bioshell_io::{read_tsv};
use bioshell_numerical::distance::euclidean_distance_squared;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Gaussian Mixture Model (GMM) estimation
/// say gmm -h to see options
struct Args {
    /// file with N-dimensional input observations: N columns of real values
    #[clap(short, long, short='i')]
    infile: String,
    /// number of clusters the points will be assigned to
    #[clap(short, long, short='k')]
    k: usize,
    /// relative error when convergence is reached
    #[clap(short, long, default_value="0.01", short='e')]
    epsilon: f64,
    /// number of times the k-means clustering is repeated; only the best solution is reported
    #[clap(short, long, default_value="1", short='r')]
    repeats: usize,
}


fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();

    let fname = args.infile;
    let sample = read_tsv(&fname);
    let n_dim: usize = sample[0].len();
    let n_data = sample.len();
    info!("{} rows loaded, data dimension is {}", n_data, n_dim);

    let k: usize = args.k;
    let mut kmeans = KMeans::new(k, sample.clone(), n_dim, euclidean_distance_squared);
    let err = kmeans.cluster_n(args.epsilon, args.repeats);
    for (v, c) in zip(&sample, kmeans.assignments()) {
        println!("{:?} {}", v, c)
    }
    info!("# Error: {err}");
}