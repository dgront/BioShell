use std::time::Instant;

use clap::{Parser};

use bioshell_cartesians::{Coordinates, CartesianSystem, coordinates_to_pdb, cubic_grid_atoms, NbList, ArgonRules, box_width};
use bioshell_sim::{Observer, ObserversSet};
use bioshell_cartesians::movers::{ChangeVolume, SingleAtomMove};
use bioshell_cartesians::observers::PdbTrajectory;
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, LennardJonesHomogenic};
use bioshell_montecarlo::{Sampler, IsothermalMC, AdaptiveMCProtocol};


#[derive(Parser, Debug)]
#[clap(name = "argon")]
#[clap(version = "0.2")]
#[clap(about = "NVT or NPT simulation of argon fluid", long_about = None)]
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// density of the system - starts NVT simulation
    #[clap(short, long, default_value_t = 0.4)]
    density: f64,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 45.0)]
    temperature: f64,
    /// pressure of an NPT simulation in [kPa]
    #[clap(short, long)] //, default_value_t = 100.0
    pressure: Option<f64>,
    /// Number of atoms in a simulation
    #[clap(short, long, default_value_t = 216)]
    natoms: usize,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
}

pub fn main() {

    // ---------- Parameters from:
    // L. A. Rowley, D. Nicholson and N. G. Parsonage
    // "Monte Carlo grand canonical ensemble calculation in a gas-liquid transition region for 12-6 Argon",
    // Journal of Computational Physics 17 pp. 401-414 (1975)
    // const EPSILON: f64 = 1.654E-21;	                // [J] per molecule
    // const EPSILON_BY_K: f64 = EPSILON / 1.381E-23; 	// = 119.8 in Kelvins
    // const SIGMA: f64 = 3.4;		                    // in Angstroms
    // ---------- Parameters from:
    // John A. White "Lennard-Jones as a model for argon and test of extended renormalization group calculations",
    // Journal of Chemical Physics 111 pp. 9352-9356 (1999)
    const EPSILON: f64 = 1.7355E-21;	                // [J] per molecule
    const EPSILON_BY_K: f64 = EPSILON / 1.380649E-23; 	// = 125.7 in Kelvins
    const SIGMA: f64 = 3.3345;		                    // in Angstroms

    const CUTOFF: f64 = 10.0;

    let args = Args::parse();
    let temperature: f64 = args.temperature;    // --- Temperature of the isothermal simulation (in the energy units)
    let n_atoms: usize = args.natoms;
    let density: f64 = args.density;
    let prefix = args.prefix;
    let tra_fname = format!("{}_tra.pdb", &prefix);
    let final_fname = format!("{}_final.pdb", &prefix);

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(SIGMA as f64,n_atoms, density));
    cubic_grid_atoms(&mut coords);
    coordinates_to_pdb(&coords,1,tra_fname.as_str(), false);

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(CUTOFF as f64,4.0,Box::new(ArgonRules{}));

    // ---------- Create the system
    let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

    // ---------- Create energy function - just the LJ term
    let lj = LennardJonesHomogenic::new(EPSILON_BY_K, SIGMA, CUTOFF);
    let pairwise_lj = PairwiseNonbondedEvaluator::new(CUTOFF as f64,lj);

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: IsothermalMC<CartesianSystem, PairwiseNonbondedEvaluator<LennardJonesHomogenic>> =
        IsothermalMC::new(temperature);
    simple_sampler.add_mover(Box::new(SingleAtomMove::new(1.0)));

    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    // ---------- If a pressure value has been provided, turn this simulation into NPT by adding a volume changing move
    if let Some(pressure) = args.pressure {
        sampler.add_mover(Box::new(ChangeVolume::new(temperature, pressure)));
    }

    // ---------- Run the simulation!
    let mut observers: ObserversSet<CartesianSystem> = ObserversSet::new();
    let mut pdb_tra = PdbTrajectory::new(tra_fname, false);
    pdb_tra.observe(&system);
    pdb_tra.if_append = true;
    observers.add_observer(Box::new(pdb_tra), 1);
    sampler.run_simulation(args.inner, args.outer, &mut system, &pairwise_lj, &mut observers);

    coordinates_to_pdb(&system.coordinates(),1,final_fname.as_str(), false);
}
