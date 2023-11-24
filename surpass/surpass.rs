use std::env;
use std::time::Instant;

use clap::{Parser};
use log::debug;
use rand::rngs::SmallRng;
use rand::SeedableRng;

use bioshell_pdb::{Structure, load_pdb_file};
use surpass::{C_N, C_O, CA_C, CaContactEnergy, CMDisplacement, ExcludedVolume, HingeMove, MoveProposal, Mover, N_CA, NonBondedEnergy, SurpassAlphaSystem, SurpassEnergy, TailMove};
use surpass::{ChainCM, RgSquared, RecordMeasurements, REndSquared};
#[allow(unused_imports)]                // NonBondedEnergyDebug can be un-commented to test non-bonded energy if the debug check fails
use surpass::{NonBondedEnergyDebug};

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
    #[clap(long, default_value = "10", short='r')]
    n_res_in_chain: usize,
    /// number of inner Monte Carlo cycles
    #[clap(long, default_value = "100", short='o')]
    outer_cycles: usize,
    /// number of inner Monte Carlo cycles
    #[clap(long, default_value = "10", short='i')]
    inner_cycles: usize,
    /// inner cycle scaling factor: the number of Monte Carlo cycles simulated for each inner MC cycle
    #[clap(long, default_value = "10", short='s')]
    cycle_factor: usize,
    /// simulation box size in Angstroms
    #[clap(long, default_value = "10000.0", short='b')]
    box_size: f64,
    /*
    /// simulation temperature
    #[clap(long, default_value = "1.0", short='t')]
    t_start: f64,
     */
}


fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // let mut strctr = load_pdb_file(&args[1]).unwrap();
    
    // --- the system
    let n_res = args.n_res_in_chain;
    if n_res < 4 { panic!("Simulated chain must be at leat 4 residues long!") }
    let n_chains = args.n_chains;
    if n_res < 4 { panic!("Simulated system must contain at least one chain!") }
    let mut rnd = SmallRng::seed_from_u64(42);
    let mut system = SurpassAlphaSystem::make_random(&vec![n_res; n_chains], args.box_size, &mut rnd);

    // --- sampling: movers and proposals
    let hinge_mover: HingeMove<4, { 4 * 4 }> = HingeMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let tail_mover_1: TailMove<1, 4> = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let mut hinge_prop: MoveProposal<4, { 4 * 4 }> = MoveProposal::new();
    let mut tail_prop_1: MoveProposal<1, 4> = MoveProposal::new();
    let tail_mover_2: TailMove<2, { 2 * 4 }> = TailMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let mut tail_prop_2: MoveProposal<2, { 2 * 4 }> = MoveProposal::new();

    // --- Observers
    let cm_measurements: Vec<ChainCM> = (0..system.count_chains()).map(|i|ChainCM::new(i)).collect();
    let mut cm =
        RecordMeasurements::new("cm.dat", cm_measurements).expect("can't write to cm.dat");
    let r2_measurements: Vec<REndSquared> = (0..system.count_chains()).map(|i|REndSquared::new(i)).collect();
    let mut rend = RecordMeasurements::new("r2.dat", r2_measurements).expect("can't write to r2.dat");
    let rg_measurements: Vec<RgSquared> = (0..system.count_chains()).map(|i|RgSquared::new(i)).collect();
    let mut rg = RecordMeasurements::new("rg.dat", rg_measurements).expect("can't write to rg.dat");
    let t_max = args.outer_cycles * args.inner_cycles / 1000;
    let mut cmd = CMDisplacement::new(system.box_length(),
                        system.count_chains(), t_max, "cm_displacement.dat");

    // --- save the starting conformation, reset the trajectory file
    // system.to_pdb_file("tra.pdb", false);

    let excl_vol = ExcludedVolume::new(&system, 3.7, 1.0);
    let energy: NonBondedEnergy<ExcludedVolume> = NonBondedEnergy::new(&system, excl_vol);
    // let mut energy: NonBondedEnergyDebug<ExcludedVolume> = NonBondedEnergyDebug::new(&system, excl_vol);

    // let cntcts = CaContactEnergy::new(&system, 10.0, -1.0, 3.7, 4.0, 5.0);
    // let energy: NonBondedEnergy<CaContactEnergy> = NonBondedEnergy::new(&system, cntcts);

    let start = Instant::now();

    let mut rnd = SmallRng::seed_from_u64(42);

    println!("initial energy: {}", energy.evaluate(&system));
    let mut _en_before: f64; // --- for debugging
    for outer in 0..args.outer_cycles {
        for inner in 0..args.inner_cycles {
            for _cycle in 0..args.cycle_factor {
                // ---------- single atom tail move
                for _tail in 0..n_chains*2 {
                    tail_mover_1.propose(&mut system, &mut rnd, &mut tail_prop_1);
                    #[cfg(debug_assertions)]
                    check_bond_lengths(&mut system, &tail_prop_1);
                    if energy.evaluate_delta(&system, &tail_prop_1) < 0.1 {
                        tail_prop_1.apply(&mut system);
                    }
                }
                // ---------- tail move of two residues
                for _tail in 0..n_chains * 2 {
                    tail_mover_2.propose(&mut system, &mut rnd, &mut tail_prop_2);
                    #[cfg(debug_assertions)]
                    check_bond_lengths(&mut system, &tail_prop_2);
                    if energy.evaluate_delta(&system, &tail_prop_2) < 0.1 {
                        tail_prop_2.apply(&mut system);
                    }
                }
                // ---------- hinge move
                for _hinge in 0..n_chains*(n_res-2) {
                    hinge_mover.propose(&mut system, &mut rnd, &mut hinge_prop);
                    let delta_e = energy.evaluate_delta(&system, &hinge_prop);
                    if delta_e < 0.1 {
                        #[cfg(debug_assertions)] {
                            _en_before = energy.evaluate(&system);
                            check_bond_lengths(&mut system, &hinge_prop);
                        }
                        hinge_prop.apply(&mut system);
                        #[cfg(debug_assertions)] {
                            let en_after = energy.evaluate(&system);
                            if (en_after - _en_before - delta_e).abs() > 0.001 {
                                // energy.report(&mut system, &hinge_prop);
                                panic!("Incorrect energy change: global {} vs delta {}\n", en_after-_en_before, delta_e);
                            }
                        }
                    }
                }
            }       // --- single inner MC cycle done (all cycle_factor MC cycles finished)
            println!("{} {} {} {:?}", outer, inner, energy.evaluate(&system), start.elapsed());
            let bond_err = system.adjust_bond_length(3.8);
            debug!("bond lengths corrected, maximum violation was: {}", bond_err);
            cm.observe(&system);
            cmd.observe(&system);
            rend.observe(&system);
            rg.observe(&system);
        }           // --- single outer MC cycle done (all inner MC cycles finished)
        // --- append a current conformation to the trajectory file
        // system.to_pdb_file("tra.pdb", true);
    }   // --- end of the simulation: all outer MC cycles done
}

#[allow(dead_code)]
fn check_bond_lengths<const N_RESIDUES: usize, const N_ATOMS: usize>(system: &mut SurpassAlphaSystem, mp: &MoveProposal<N_RESIDUES, N_ATOMS>) {
    let mut backup: MoveProposal<N_RESIDUES, N_ATOMS> = MoveProposal::new();
    backup.first_moved_pos = mp.first_moved_pos;
    backup.backup(system);
    mp.apply(system);
    let true_bond_lengths = [C_N, N_CA, CA_C, C_O];
    for i_residue in 0..system.count_residues()-1 {
        let mut idx = i_residue*4;
        // --- check C-N bond only when both atoms are of the same chain; also, the first residue doesn't have preceding C
        if system.chain(i_residue) == system.chain(i_residue +1) && i_residue > 0 {
            check_single_bond(system, idx-1, idx, true_bond_lengths[0], &backup);
        }
        for i_atom in 1..4 {
            idx += 1;
            check_single_bond(system, idx-1, idx, true_bond_lengths[i_atom], &backup);
        }
    }
    backup.apply(system);
}

fn check_single_bond<const N_RESIDUES: usize, const N_ATOMS: usize>(system: &mut SurpassAlphaSystem,
            i_atom: usize, j_atom: usize, d_expected: f64, backup: &MoveProposal<N_RESIDUES, N_ATOMS>) {
    let d_observed = system.distance(i_atom, j_atom);
    let i_residue = i_atom / 4;
    if (d_observed - d_expected).abs() > 0.01 {
        system.to_pdb_file("after.pdb", false);
        let av = system.atom_to_vec3(i_residue + 1);
        let prev_av = system.atom_to_vec3(i_residue);
        backup.apply(system);
        let bv = system.atom_to_vec3(i_residue + 1);
        let prev_bv = system.atom_to_vec3(i_residue);
        system.to_pdb_file("before.pdb", false);

        panic!("Broken bond between {} and {} of {}, current length is: {}\nPos. before: {} {}\nPos. after: {} {}\n",
               i_atom, j_atom, i_residue, d_observed, &prev_bv, &bv, &prev_av, &av);
    }
}

#[allow(dead_code)]
fn check_delta_en(en_before: f64, en_after: f64, delta: f64) {
    if (en_after-en_before-delta).abs() > 0.001 {
        panic!("Incorrect energy change: global {} vs delta {}\n", en_after-en_before, delta);
    }
}