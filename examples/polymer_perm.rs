use std::{env, fmt};
use std::time::Instant;

use clap::{Parser};
use bioshell_cartesians::{CartesianSystem, Coordinates, NbList, PolymerRules};
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};

use bioshell_sim::{Energy, System};


#[derive(Parser, Debug)]
#[clap(name = "ising")]
#[clap(about = "PERM generator for simple polymer chains", long_about = None)]
struct Args {
    /// loads a staring conformation from a text file
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.70)]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 100, short='n')]
    n_beads: usize,
    /// number of chains to be generated for each cycle
    #[clap(short='k', long, default_value_t = 100)]
    chains: usize,
    /// number of PERM cycles to be performed; at every cycle n_chains are created and W+ and W- bounding parameters are updated
    #[clap(short='c', long, default_value_t = 100)]
    cycles: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
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
    let temperature: f64 = args.temperature;    // --- Temperature of the isothermal simulation (in the energy units)
    let prefix = args.prefix;
    let tra_fname = format!("{}_tra.txt", &prefix);

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(args.n_beads);
    coords.set_box_len(1000.0);
    let nbl = NbList::new(E_TO, 1.0, Box::new(PolymerRules{}));
    let system = CartesianSystem::new(coords, nbl);

    // ---------- Create an energy function - just contacts for now
    // ---------- Contact energy
    let contact_kernel = SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL);
    let contacts: PairwiseNonbondedEvaluator<SimpleContact> = PairwiseNonbondedEvaluator::new(E_TO as f64, contact_kernel);

    // ---------- Create a PERM generator

    // ---------- Run the generator!
    let start = Instant::now();
    for i_cycle in 0..args.cycles {
        for i_chain in 0..args.chains {
        }
    }
}
