use clap::{Parser};
use bioshell_numerical::clustering::{EuclideanNeighbors, Optics};
use bioshell_core::utils::{read_tsv, out_writer};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Gaussian Mixture Model (GMM) estimation
/// say gmm -h to see available options
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


fn main() {

    let args = Args::parse();

    // ---------- input data
    let sample = read_tsv(&args.infile);
    println!("{} rows loaded, data dimension is {}", sample[0].len(), sample.len());

    // ---------- clustering parameters
    let min_points: usize = args.min_points;
    let min_cluster: usize = args.min_cluster;
    let epsilon: f64 = args.epsilon;

    let opt_clust = Optics::new(epsilon, min_points,
                                    Box::new(EuclideanNeighbors::new(sample.clone())));

    let mut clusters = opt_clust.clusters();
    clusters.sort_by(|c1, c2| c2.len().partial_cmp(&c1.len()).unwrap());
    let mut i:i16 = 0;
    for ci in &clusters {
        if ci.len() < min_cluster { continue; }
        i += 1;
        let mut out = out_writer(format!("c{}.dat", i).as_str());
        for iel in ci {
            for val in &sample[*iel] { write!(out, "{}\t", val); }
            write!(out, "\n");
        }
    }

}