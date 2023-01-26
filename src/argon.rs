use std::time::Instant;

use clap::{Parser};

use bioshell_cartesians::{Coordinates, CartesianSystem, coordinates_to_pdb, cubic_grid_atoms, NbList, ArgonRules, VolumeChangingProtocol};
use bioshell_sim::{Energy, System};
use bioshell_cartesians::movers::SingleAtomMove;
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, LennardJonesHomogenic};
use bioshell_montecarlo::{Sampler, AcceptanceStatistics, MCProtocol, MetropolisCriterion,
                          AdaptiveMCProtocol, MoversSet};


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

fn box_width(sphere_radius: f64, n_spheres: usize, density: f64) -> f64 {

    let v: f64 = 4.0/3.0 * std::f64::consts::PI * sphere_radius.powi(3);
    (n_spheres as f64 * v/density).powf(1.0/3.0)
}

pub fn main() {

    const EPSILON: f64 = 1.654E-21;	                // [J] per molecule
    const EPSILON_BY_K: f64 = EPSILON / 1.381E-23; 	// = 119.6 in Kelvins
    const SIGMA: f64 = 3.4;		                    // in Angstroms
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
    // let energy: Box<dyn Energy<CartesianSystem>> = Box::new(pairwise);

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: MCProtocol<CartesianSystem, PairwiseNonbondedEvaluator<LennardJonesHomogenic>, MetropolisCriterion> =
        MCProtocol::new(MetropolisCriterion::new(temperature));
    simple_sampler.add_mover(Box::new(SingleAtomMove::new(1.0)));
    let mut pressure = 1.0;
    if let Some(p) = args.pressure {
        pressure = p;
    }
    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    // let mut sampler = VolumeChangingProtocol::new(pressure, Box::new(sampler));

    // ---------- Run the simulation!
    let start = Instant::now();
    let mut recent_acceptance = AcceptanceStatistics::default();
    for i in 0..args.outer {
        let stats = sampler.get_mover(0).acceptance_statistics();
        sampler.make_sweeps(args.inner,&mut system, &pairwise_lj);
        coordinates_to_pdb(&system.coordinates(), (i+1) as i16, tra_fname.as_str(), true);
        let f_succ = stats.recent_success_rate(&recent_acceptance);
        recent_acceptance = stats;
        println!("{:6} {:9.3} {:5.3} {:.2?}", i, pairwise_lj.energy(&system) / system.size() as f64, f_succ, start.elapsed());
    }
    coordinates_to_pdb(&system.coordinates(),1,final_fname.as_str(), false);
}
