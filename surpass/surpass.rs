use std::env;

use clap::{Parser};

use bioshell_pdb::{Structure, load_pdb_file};
use surpass::{ExcludedVolume, HingeMove, MoveProposal, Mover, SurpassAlphaSystem, SurpassEnergy, TailMove};

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

    let excl_vol = ExcludedVolume::new(&system, 3.7, 1.0);
    println!("{}", excl_vol.evaluate(&system));
    for outer in 0..args.outer_cycles {
        for inner in 0..args.inner_cycles {
            for _tail in 0..n_chains*2 {
                tail_mover.propose(&mut system, &mut tail_prop);
                #[cfg(debug_assertions)]
                check_bond_lengths(&mut system, &tail_prop, 3.8);
                if excl_vol.evaluate_delta(&system, &tail_prop) < 0.1 {
                    tail_prop.apply(&mut system);
                }
            }
            for _hinge in 0..n_chains*(n_res-2) {
                let en_before = excl_vol.evaluate(&system);
                hinge_mover.propose(&mut system, &mut hinge_prop);
                #[cfg(debug_assertions)]
                check_bond_lengths(&mut system, &hinge_prop, 3.8);
                let delta_e = excl_vol.evaluate_delta(&system, &hinge_prop);
                if delta_e < 0.1 {
                    hinge_prop.apply(&mut system);
                    let en_after = excl_vol.evaluate(&system);
                    #[cfg(debug_assertions)]
                    check_delta_en(en_before, en_after, delta_e);
                }
            }
        }
        println!("{} {}", outer, excl_vol.evaluate(&system));
        // --- append a current conformation to the trajectory file
        system.to_pdb_file("tra.pdb", true);
    }
}

fn check_bond_lengths<const N: usize>(system: &mut SurpassAlphaSystem, mp: &MoveProposal<N>, d: f64) {
    let mut backup: MoveProposal<N> = MoveProposal::new();
    backup.first_moved_pos = mp.first_moved_pos;
    backup.backup(system);
    mp.apply(system);
    for i in 0..system.count_atoms()-1 {
        if system.chain(i) != system.chain(i+1) { continue }
        let dd = system.distance(i+1, i);
        if (dd-d).abs() > 0.01 {
            system.to_pdb_file("after.pdb", false);
            let av = system.ca_to_vec3(i+1);
            let prev_av = system.ca_to_vec3(i);
            backup.apply(system);
            let bv = system.ca_to_vec3(i+1);
            let prev_bv = system.ca_to_vec3(i);
            system.to_pdb_file("before.pdb", false);

            panic!("Broken bond between {} and {}, current length is: {}\nPos. before: {} {}\nPos. after: {} {}\n",
                   i, i+1, dd, &prev_bv, &bv, &prev_av, &av);
        }
    }
    backup.apply(system);
}

fn check_delta_en(en_before: f64, en_after: f64, delta: f64) {
    if (en_after-en_before-delta).abs() > 0.001 {
        panic!("Incorrect energy change: global {} vs delta {}\n", en_after-en_before, delta);
    }
}