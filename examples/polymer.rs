use std::env;
use std::time::Instant;

use clap::{Parser};
use log::{info};

use bioshell_cartesians::{Coordinates, CartesianSystem, coordinates_to_pdb, RandomChain, pdb_to_coordinates,
                          NbList, PolymerRules, box_width};
use bioshell_cartesians::movers::SingleAtomMove;
use bioshell_cartesians::observers::{REndSquared, GyrationSquared, PdbTrajectory};

use bioshell_sim::{Energy, ObserversSet, System};

use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{Sampler, IsothermalMC, AdaptiveMCProtocol, StepwiseBuilder};

use bioshell_ff::{TotalEnergy};
use bioshell_ff::bonded::SimpleHarmonic;

#[derive(Parser, Debug)]
#[clap(name = "polymer")]
#[clap(version = "0.2")]
#[clap(about = "Simple polymer model", long_about = None)]
/// simulation of a simple polymer model
/// say polymer -h to see options
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.0, short='t')]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 100, short='n')]
    nbeads: usize,
    /// by-volume density of the simulated system
    #[clap(short, long, default_value_t = 0.4, short='d')]
    density: f64,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
    /// thickness of the buffer zone used by non-bonded list to hash non-bonded interactions
    #[clap(long, default_value_t = 4.0)]
    buffer: f64,
}


pub fn main() {
    const E_REP: f64 = 4.25;
    const E_FROM: f64 = 4.5;
    const E_TO: f64 = 6.0;
    const E_VAL: f64 = -1.0;
    const REP_VAL: f64 = 1000.0;
    const MAX_MOVE_RANGE: f64 = 1.0;

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // ---------- Temperature of the isothermal simulation (in the energy units)
    let temperature: f64 = args.temperature;    // --- Temperature of the isothermal simulation (in the energy units)
    let density: f64 = args.density;            // --- density of the system
    let box_length: f64;
    let prefix = args.prefix;                   // --- prefix for the output files
    let tra_fname = format!("{}tra.pdb", &prefix);

    // ---------- NBL buffer radius
    let buffer_thickness = args.buffer.max(MAX_MOVE_RANGE * 4.0);

    // ---------- Create system's coordinates
    let n_beads: usize;
    let mut coords: Coordinates;

    if args.infile == "" {
        n_beads = args.nbeads;
        coords = Coordinates::new(n_beads);
        box_length = box_width(E_REP, n_beads, density);
        coords.set_box_len(box_length);
    } else {
        let res = pdb_to_coordinates(&args.infile).ok();
        match res {
            Some(x) => { coords = x; },
            _ => panic!("Can't read from file >{}<", args.infile)
        };
        n_beads = coords.size();
        box_length = box_width(E_REP, n_beads, density);
        coords.set_box_len(box_length);
    }

    info!("Simulation settings:
        temperature: {:.3}
        no. beads:   {}
        box length:  {:.3}
        NBL buffer:  {:.3}", temperature, n_beads, box_length, buffer_thickness);

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(E_TO,buffer_thickness,Box::new(PolymerRules{}));
    // ---------- Create the system
    let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

    // ---------- Contact energy
    let contact_kernel = SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL);
    let contacts: PairwiseNonbondedEvaluator<SimpleContact> = PairwiseNonbondedEvaluator::new(E_TO as f64, contact_kernel);
    // ---------- Harmonic energy (i.e. springs between beads)
    let harmonic = SimpleHarmonic::new(3.8,1.0);
    // ---------- Total energy contains the contacts energy and bond springs
    let mut total = TotalEnergy::new();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    if system.size() == 0 { RandomChain::default().build(&mut system, &total); }

    // ---------- Show the starting point
    coordinates_to_pdb(&system.coordinates(), 0, tra_fname.as_str(), false);
    println!("# starting energy: {}", total.energy(&system));

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: IsothermalMC<CartesianSystem, TotalEnergy<CartesianSystem>> = IsothermalMC::new(temperature);
    simple_sampler.add_mover(Box::new(SingleAtomMove::new(MAX_MOVE_RANGE)));

    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    let mut observations: ObserversSet<CartesianSystem> = ObserversSet::new();
    observations.add_observer(Box::new(PdbTrajectory::new(tra_fname, true)), 1);
    observations.add_observer(Box::new(GyrationSquared::new("rg.dat".to_string(), true)), 1);
    observations.add_observer(Box::new(REndSquared::new("r2.dat".to_string(), true)), 1);
    sampler.run_simulation(args.inner, args.outer, &mut system, &total, &mut observations);

    coordinates_to_pdb(&system.coordinates(), 1i16, format!("{}final.pdb", &prefix).as_str(), false);
}
