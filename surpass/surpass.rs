use std::env;

use clap::{Parser};

use bioshell_pdb::{Structure, load_pdb_file};
use surpass::{CaContactEnergy, ExcludedVolume, HingeMove, MoveProposal, Mover, NonBondedEnergy, SurpassAlphaSystem, SurpassEnergy, TailMove};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Surpass-alpha model for coarse grained simulations of proteins
/// say surpass -h to see options
struct Args {
    /// starting conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// number of chains of the simulated system
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
    let energy: NonBondedEnergy<ExcludedVolume> = NonBondedEnergy::new(&system, excl_vol.repulsion_cutoff(), excl_vol);

    // let cntcts = CaContactEnergy::new(&system, 10.0, -1.0, 3.7, 4.0, 5.0);
    // let energy: NonBondedEnergy<CaContactEnergy> = NonBondedEnergy::new(&system, 3.7, cntcts);


    println!("initial energy: {}", energy.evaluate(&system));
    let mut _en_before: f64; // --- for debugging
    for outer in 0..args.outer_cycles {
        for inner in 0..args.inner_cycles {
            for _tail in 0..n_chains*2 {
                tail_mover.propose(&mut system, &mut tail_prop);
                #[cfg(debug_assertions)]
                check_bond_lengths(&mut system, &tail_prop, 3.8);
                if energy.evaluate_delta(&system, &tail_prop) < 0.1 {
                    tail_prop.apply(&mut system);
                }
            }
            for _hinge in 0..n_chains*(n_res-2) {
                hinge_mover.propose(&mut system, &mut hinge_prop);
                let delta_e = energy.evaluate_delta(&system, &hinge_prop);
                if delta_e < 0.1 {
                    #[cfg(debug_assertions)] {
                        _en_before = energy.evaluate(&system);
                        check_bond_lengths(&mut system, &hinge_prop, 3.8);
                    }
                    hinge_prop.apply(&mut system);
                    #[cfg(debug_assertions)] {
                        let en_after = energy.evaluate(&system);
                        check_delta_en(_en_before, en_after, delta_e);
                    }
                }
            }
        }
        println!("{} {}", outer, energy.evaluate(&system));
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