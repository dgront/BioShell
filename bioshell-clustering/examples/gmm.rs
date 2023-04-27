use std::env;
use std::iter::zip;
use nalgebra::{DMatrix, DVector};
use clap::{Parser};
use log::{info};

use bioshell_statistics::{MultiNormalDistribution, OnlineMultivariateStatistics};
use bioshell_clustering::em::expectation_maximization;
use bioshell_clustering::kmeans::KMeans;
use bioshell_core::utils::{read_tsv};
use bioshell_numerical::distance::euclidean_distance_squared;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Gaussian Mixture Model (GMM) estimation
/// say gmm -h to see options
struct Args {
    /// file with N-dimensional input observations: N columns of real values
    #[clap(short, long, default_value = "", short='i')]
    infile: String,
    /// number of N-dimensional normal distributions to be inferred
    #[clap(short, long, default_value = "3", short='k')]
    n_distr: usize,
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

    let n_dist: usize = args.n_distr;

    let mut kmeans = KMeans::new(n_dist, sample.clone(), n_dim, euclidean_distance_squared);
    let err = kmeans.cluster_n(0.001, 100);
    info!("k-means error: {}", err);
    info!("k-means cluster sizes: {:?}", &kmeans.sizes());

    let mut stats: Vec<OnlineMultivariateStatistics> = vec![OnlineMultivariateStatistics::new(n_dim); n_dist];
    for (v, c) in zip(&sample, kmeans.assignments()) {
        stats[*c].accumulate(v);
    }

    // ---------- Initial multivariate normal distributions based on K-means results
    let mut normals = vec!(MultiNormalDistribution::new(n_dim); n_dist);
    for i in 0..n_dist {
        let mean = stats[i].avg();
        let mut cov_val: Vec<f64> = vec![];
        for v in stats[i].cov() {
            for vi in v { cov_val.push(vi); }
        }
        normals[i].set_parameters(&DVector::from_vec(mean.clone()),
            &DMatrix::from_vec(n_dim, n_dim, cov_val));
    }

    let mut weights = vec![0.0;n_dist];
    let log_likelihood = expectation_maximization(&mut normals, &sample, &mut weights, 0.000001);
    info!("Log-likelihood average: {}", log_likelihood / n_data as f64);
    println!("[");
    for i in 0..n_dist {
        println!("\t{{'n': {:5}, {} }},", weights[i], &normals[i]);
    }
    println!("]");
}