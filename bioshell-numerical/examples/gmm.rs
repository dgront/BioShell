use std::fs::File;
use nalgebra::{DMatrix, DVector};
use clap::{Parser};
use csv;
use rand::Rng;

use bioshell_numerical::statistics::{MultiNormalDistribution, expectation_maximization};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Gaussian Mixture Model (GMM) estimation
/// say gmm -h to see options
struct Args {
    /// input observations
    #[clap(short, long, default_value = "", short='i')]
    infile: String,
    #[clap(short, long, default_value = "3", short='n')]
    n_distr: usize,
}

fn read_csv(fname:&str) -> Vec<Vec<f64>> {

    let file = File::open(fname).unwrap();
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(file);
    let mut data : Vec<Vec<f64>> = Vec::new();
    for record in rdr.records() {
        if let Ok(r) = record {
            let mut row : Vec<f64> = r.iter().map(|e| e.parse::<f64>().unwrap()).collect();
            data.push(row);
        }
    }

    return data;
}

fn main() {

    let args = Args::parse();

    let fname = args.infile;
    let sample = read_csv(&fname);
    let n_dim: usize = sample[0].len();
    let n_data = sample.len();
    println!("{} rows loaded, data dimension is {}", n_data, n_dim);

    let n_dist: usize = args.n_distr;


    // ---------- Distributions to be inferred
    let mut normals = vec!(MultiNormalDistribution::new(n_dim); n_dist);
    // ---------- starting parameters
    let mean: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 0.0, 0.0],
                                   vec![1.35, 1.8, 0.0, 0.0],
                                   vec![1.35, 2.5, 0.0, 0.0]];
    let sdev_diag: Vec<Vec<f64>> = vec![vec![0.2, 0.2, 1.0, 1.0],
                                   vec![0.2, 0.2, 1.0, 1.0],
                                   vec![0.2, 0.2, 1.0, 1.0]];

    for i in 0..n_dist {
        normals[i].set_parameters(&DVector::from_vec(mean[i].clone()),
            &DMatrix::from_diagonal(&DVector::from_vec(sdev_diag[i].clone())));
    }

    // ---------- Randomly assign data points to distributions
    let mut assignment: Vec<Vec<bool>> = vec!(vec!(false; n_data); n_dist);
    for i in 0..n_data {
        let i_dist = rand::thread_rng().gen_range(0..n_dist);
        assignment[i_dist][i] = true;
    }

    expectation_maximization(&mut normals, &sample, &mut assignment, 0.000001);
    let mut totals: Vec<i32> = vec!(0; n_dist);
    for i in 0..n_dist {
        totals[i] = (&assignment[i]).into_iter().map(|v| *v as i32).sum();
        println!("{} {}",totals[i], &normals[i]);
    }
}