use std::env;
use std::io::Error;
use clap::{Parser};
use bioshell_clustering::{CartesianPoints};
use bioshell_clustering::optics::{Optics};
use bioshell_io::{read_delimited_values, out_writer, open_file};

use log::{error, info};
use bioshell_datastructures::euclidean_distance_squared;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
///  OPTICS (Ordering Points To Identify the Clustering Structure) clustering algorithm
/// say optics -h to see available options
struct Args {
    /// input observations
    #[clap(short, long, short='i')]
    infile: String,
    #[clap(short, long, default_value = "5", short='n')]
    min_points: usize,
    #[clap(short, long, default_value = "2", short='k')]
    min_cluster: usize,
    #[clap(short, long, short='e')]
    epsilon: f64,
}

/// A simple application that performs clustering with the OPTICS method
/// USAGE:
///     optics -h
fn main() -> Result<(), Error> {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();

    // ---------- input data
    let fname = &args.infile;
    let reader = open_file(fname).expect(&format!("Can't open {} file!", fname));
    let sample = read_delimited_values(reader, b'\t')?;
    info!("{} rows loaded, data dimension is {}", sample.len(), sample[0].len());

    // ---------- clustering parameters
    let min_points: usize = args.min_points;
    let min_cluster: usize = args.min_cluster;
    let epsilon: f64 = args.epsilon;

    let opt_clust = Optics::new(epsilon, min_points,
        Box::new(CartesianPoints::new(euclidean_distance_squared, sample.clone(), sample[0].len())));

    let mut clusters = opt_clust.clusters();
    clusters.sort_by(|c1, c2| c2.len().partial_cmp(&c1.len()).unwrap());
    let mut i:i16 = 0;
    for ci in &clusters {
        if ci.len() < min_cluster { continue; }
        i += 1;
        let fname = format!("c{}.dat", i);
        let mut out = out_writer(fname.as_str(), false);
        for iel in ci {
            match write!(out, "{}", sample[*iel][0]) {
                Err(e) => error!("Error while writing to {} file: {:?}", fname, e),
                _ => ()
            }
            for ival in 1..sample[*iel].len() {
                match write!(out, "\t{}", sample[*iel][ival]) {
                    Err(e) => error!("Error while writing to {} file: {:?}", fname, e),
                    _ => ()
                }
            }

            write!(out, "\n").ok();
        }
    }

    Ok(())
}