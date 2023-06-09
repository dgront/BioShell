use clap::Parser;
use log::{error, info};
use rand::Rng;
use std::env;
use std::time::Instant;

use bioshell_cartesians::{coordinates_to_pdb, Coordinates};
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{StepwiseBuilder, StepwiseMover, PERM};
use bioshell_numerical::Vec3;

use bioshell_sim::{Energy, ResizableSystem, System};

/// Builds a simple chain on 2D square lattice
pub struct ChainBuilder {
    /// temperature for the simulation
    pub temperature: f64,
    /// bond length - use 4.0 for nicely looking PDB trajectory
    pub bond_length: f64,
    weights: [f64; 4],      // scratch memory to store weight for a single step
    moves: [(f64, f64); 4], // scratch memory to store proposed move coordinates
}

impl ChainBuilder {
    /// Create a new chain builder that construct chains from Boltzmann distribution at given temperature
    pub fn new(temperature: f64) -> ChainBuilder {
        ChainBuilder {
            temperature,
            bond_length: 4.0,
            weights: [0.0, 0.0, 0.0, 0.0],
            moves: [(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)],
        }
    }
}

impl<E: Energy<Coordinates>> StepwiseMover<Coordinates, E> for ChainBuilder {
    /// Starts a new chain by placing its first bead in (0, 0, 0) and returns 1.0 as the statistical weight
    fn start(&mut self, system: &mut Coordinates, _energy: &E) -> f64 {
        let c = system.get_box_len() / 2.0;
        system.set_xyz(0, c, c, 0.0);
        system.set_size(1);

        return 1.0;
    }

    /// Grows the current chain by one bead
    fn grow_by_one(&mut self, system: &mut Coordinates, energy: &E) -> f64 {
        let i = system.get_size();
        system.set_size(i + 1);

        let center: Vec3 = (&system[i - 1]).clone();
        // ---------- propose n_trials random proposals and score them
        for k in 0..4 {
            system.set_xyz(
                i,
                center.x + self.moves[k].0 * self.bond_length,
                center.y + self.moves[k].1 * self.bond_length,
                0.0,
            );
            let en = energy.energy_by_pos(system, i);
            self.weights[k] = (-en / self.temperature).exp();
        }
        let total = self.weights.iter().sum();
        if total < 1e-100 {
            return 0.0;
        } // --- no suitable move generated

        // ---------- select one of the possible extension by importance sampling
        let mut rng = rand::thread_rng();
        let r = rng.gen_range(0.0..total);
        let mut which_v: usize = 0;
        let mut s = self.weights[which_v];
        while s <= r {
            which_v += 1;
            s += self.weights[which_v]
        }

        // ---------- set the coordinates
        system.set_xyz(
            i,
            center.x + self.moves[which_v].0 * self.bond_length,
            center.y + self.moves[which_v].1 * self.bond_length,
            0.0,
        );

        // ---------- return the statistical weight
        return total;
    }
}

#[derive(Parser, Debug)]
#[clap(name = "chain2d_perm")]
#[clap(about = "Simple PERM demo builds chains on a 2D square lattice", long_about = None)]
struct Args {
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.0)]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 25, short = 'n')]
    n_beads: usize,
    /// number of chains to be generated
    #[clap(short = 'k', long, default_value_t = 10)]
    chains: usize,
}

pub fn main() {
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Args::parse();
    let temperature: f64 = args.temperature; // --- Temperature of the isothermal simulation (in the energy units)

    // ---------- Create system's coordinates
    let mut system = Coordinates::new(args.n_beads);
    system.set_box_len(1000.0);

    // ---------- Create an energy function - just contacts
    let bond_length = 4.0;
    let contact_kernel =
        SimpleContact::new(0.1, 0.9 * bond_length, 1.0 * bond_length, 10000.0, -1.0);
    let contacts: PairwiseNonbondedEvaluator<SimpleContact> =
        PairwiseNonbondedEvaluator::new(1.1 * bond_length, contact_kernel);

    // ---------- Create a PERM generator
    let mover = ChainBuilder::new(temperature);
    let mut sampler = PERM::new(args.n_beads, 0.01, 1.0, Box::new(mover));
    // ---------- Run the generator!
    let start = Instant::now();
    let mut i_chain = 0;
    while i_chain < args.chains {
        let w = sampler.build(&mut system, &contacts);
        if w > 0.0 {
            i_chain += 1;
            println!("en, w: {:5} {:5.1e}", contacts.energy(&system), w);
            coordinates_to_pdb(&system, i_chain as i16, "chain2d_tra.pdb", true);
            info!("chain {:5} finished after {:.2?}", i_chain, start.elapsed());
            // println!("{}", sampler);
        }
    }
}
