use std::env;

use clap::{Parser};

use bioshell_pdb::{Structure, load_pdb_file};
use surpass::{HingeMove, MoveProposal, Mover, SurpassAlphaSystem, TailMove};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Gaussian Mixture Model (GMM) estimation
/// say gmm -h to see options
struct Args {
    /// file with N-dimensional input observations: N columns of real values
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// number of N-dimensional normal distributions to be inferred
    #[clap(short, long, default_value = "1", short='c')]
    n_chains: usize,
    /// number of atoms in every chain of the system; the total number of residues is n_chains x n_res_in_chain
    #[clap(long, default_value = "100", short='r')]
    n_res_in_chain: usize,
    /// number of inner Monte Carlo cycles
    #[clap(long, default_value = "100", short='i')]
    inner_cycles: usize,
    /// number of inner Monte Carlo cycles
    #[clap(long, default_value = "100", short='o')]
    outer_cycles: usize,
    /// simulation box size in Angstroms
    #[clap(long, default_value = "10000.0", short='b')]
    box_size: f64,
    /// simulation temperature
    #[clap(long, default_value = "1.0", short='t')]
    t_start: f64,
}


fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // let mut strctr = load_pdb_file(&args[1]).unwrap();

    let n_res = args.n_res_in_chain;
    let n_chains = args.n_chains;
    let mut system = SurpassAlphaSystem::new(&vec![n_res; n_chains], args.box_size);

    let hinge_mover: HingeMove<4> = HingeMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let tail_mover: TailMove = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let mut hinge_prop: MoveProposal<4> = MoveProposal::new();
    let mut tail_prop: MoveProposal<1> = MoveProposal::new();

    // --- save the starting conformation, reset the trajectory file
    system.to_pdb_file("tra.pdb", false);

    for _outer in 0..args.outer_cycles {
        for _inner in 0..args.inner_cycles {
            for _tail in 0..n_chains*2 {
                tail_mover.propose(&mut system, &mut tail_prop);
                tail_prop.apply(&mut system);
            }
            for _hinge in 0..n_chains*(n_res-2) {
                hinge_mover.propose(&mut system, &mut hinge_prop);
                hinge_prop.apply(&mut system);
            }
        }
        // --- append a current conformation to the trajectory file
        system.to_pdb_file("tra.pdb", true);
    }
}