use std::env;
use std::io::Error;
use nalgebra::{DMatrix, DVector};
use clap::{Parser};
use log::{debug, info};

use bioshell_statistics::{MultiNormalDistribution, OnlineMultivariateStatistics};
use bioshell_clustering::em::expectation_maximization;
use bioshell_clustering::kmeans::KMeans;
use bioshell_datastructures::euclidean_distance_squared;
use bioshell_io::{open_file, read_delimited_values};


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
    /// number of K-means repeats to get initial data partition
    #[clap(long, default_value = "10")]
    kmeans_repeats: usize,
    /// number of times the Gaussian mixture is fit to the data; the best solution will be reported
    #[clap(long, default_value = "100")]
    gmm_repeats: usize,
}

fn run_once(args: &Args, sample: &mut Vec<Vec<f64>>,
            gmm: &mut Vec<MultiNormalDistribution>, weights: &mut Vec<f64>) -> f64 {

    assert_eq!(gmm.len(), weights.len());

    let n_dim: usize = sample[0].len();
    let n_dist = gmm.len();

    let mut kmeans = KMeans::new(n_dist, sample.clone(), n_dim, euclidean_distance_squared);
    let err = kmeans.cluster_n(0.001, args.kmeans_repeats);
    info!("k-means lowest error: {}", err);
    info!("k-means cluster sizes: {:?}", &kmeans.sizes());

    // ---------- Initial multivariate normal distributions based on K-means results
    let mut stats: Vec<OnlineMultivariateStatistics> = vec![OnlineMultivariateStatistics::new(n_dim); n_dist];

    let point_assignment = kmeans.assignments();
    for i in 0..sample.len() {
        stats[point_assignment[i]].accumulate(&sample[i]);
    }
    // for (v, c) in zip(sample, kmeans.assignments()) {
    //     stats[*c].accumulate(v);
    // }
    for i in 0..n_dist {
        let mean = stats[i].avg();
        let mut cov_val: Vec<f64> = vec![];
        for v in stats[i].cov() {
            for vi in v { cov_val.push(vi); }
        }
        gmm[i].set_parameters(&DVector::from_vec(mean.clone()),
                                  &DMatrix::from_vec(n_dim, n_dim, cov_val));
    }

    let likelihood = expectation_maximization(gmm, sample, weights, 0.000001);

    return likelihood;
}

fn main() -> Result<(), Error> {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();

    let fname = &args.infile;
    let reader = open_file(fname).expect(&format!("Can't open {} file!", fname));
    let mut sample = read_delimited_values(reader, b'\t')?;
    let n_dim: usize = sample[0].len();
    let n_data = sample.len();
    info!("{} rows loaded, data dimension is {}", n_data, n_dim);

    let n_dist: usize = args.n_distr;
    let mut normals = vec!(MultiNormalDistribution::new(n_dim); n_dist);
    let mut weights = vec![0.0; n_dist];

    let mut best_likelihood = run_once(&args, &mut sample,&mut normals, &mut weights);
    info!("New best log-likelihood average: {}", best_likelihood / n_data as f64);
    let mut best_normals = normals.clone();
    let mut best_weights = weights.clone();
    for _k in 1..args.gmm_repeats {
        let log_likelihood = run_once(&args, &mut sample,&mut normals, &mut weights);
        if best_likelihood < log_likelihood {
            best_normals = normals.clone();
            best_weights = weights.clone();
            best_likelihood = log_likelihood;
            info!("New best log-likelihood average: {}", best_likelihood / n_data as f64);
        }
    }
    info!("Best log-likelihood average: {}", best_likelihood / n_data as f64);
    println!("[");
    for i in 0..n_dist {
        println!("\t{{'n': {:5}, {} }},", best_weights[i], &best_normals[i]);
    }
    println!("]");

    Ok(())
}