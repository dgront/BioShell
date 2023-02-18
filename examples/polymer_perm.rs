use std::{env};
use std::time::Instant;
use log::{info, error};

use clap::{Parser};
use bioshell_cartesians::{CartesianSystem, Coordinates, coordinates_to_pdb, PERMChainStep};
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{PERM, StepwiseBuilder};

use bioshell_sim::{Energy};

#[derive(Parser, Debug)]
#[clap(name = "polymer_perm")]
#[clap(about = "PERM generator for simple polymer chains", long_about = None)]
struct Args {
    /// loads a staring conformation from a text file
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.70)]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 20, short='n')]
    n_beads: usize,
    /// number of chains to be generated
    #[clap(short='k', long, default_value_t = 10)]
    chains: usize,
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

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    let temperature: f64 = args.temperature;    // --- Temperature of the isothermal simulation (in the energy units)
    let prefix = args.prefix;
    let tra_fname = format!("{}_tra.pdb", &prefix);

    // ---------- Create system's coordinates
    let mut system = Coordinates::new(args.n_beads);
    system.set_box_len(1000.0);

    // ---------- Create an energy function - just contacts for now
    // ---------- Contact energy
    let contact_kernel = SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL);
    let contacts: PairwiseNonbondedEvaluator<SimpleContact> = PairwiseNonbondedEvaluator::new(E_TO as f64, contact_kernel);

    // ---------- Create a PERM generator
    let mover = PERMChainStep::new(temperature, 16);
    let mut sampler = PERM::new(args.n_beads,0.1, 10.0,  Box::new(mover));
    // ---------- Run the generator!
    let start = Instant::now();
    let mut i_chain = 0;
    while i_chain < args.chains {
        let w = sampler.build(&mut system, &contacts);
        if w > 0.0 {
            i_chain += 1;
            println!("en, w: {:5} {:5.1e}", contacts.energy(&system), w);
            // coordinates_to_pdb(&system, i_chain as i16, tra_fname.as_str(), true);
        }

        info!("Chain {:7} finished after {:.2?}", i_chain, start.elapsed());
        // println!("{}", sampler);
    }
}
