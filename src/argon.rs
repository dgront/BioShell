use std::time::Instant;
use std::ops::Range;

use clap::Parser;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb};
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, LennardJonesHomogenic};
use bioshell_sim::generators::cubic_grid_atoms;
use bioshell_sim::sampling::movers::{single_atom_move};
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// density of the system
    #[clap(short, long, default_value_t = 0.4)]
    density: f32,
    /// density of the system
    #[clap(short, long, default_value_t = 45.0)]
    temperature: f32,
    /// Number of atoms in a simulation
    #[clap(short, long, default_value_t = 216)]
    natoms: usize,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
}

fn box_width(sphere_radius: f32, n_spheres: usize, density: f32) -> f32 {

    let v: f32 = 4.0/3.0 * std::f32::consts::PI * sphere_radius.powi(3);
    (n_spheres as f32 * v/density).powf(1.0/3.0)
}

pub fn main() {

    const EPSILON: f64 = 1.654E-21;	                // [J] per molecule
    const EPSILON_BY_K: f64 = EPSILON / 1.381E-23; 	// = 119.6 in Kelvins
    const SIGMA: f64 = 3.4;		                    // in Angstroms
    const CUTOFF: f64 = 10.0;

    let args = Args::parse();

    let n_atoms: usize = args.natoms;
    let density: f32 = args.density;
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(SIGMA as f32,n_atoms, density));
    cubic_grid_atoms(&mut coords);
    to_pdb(&coords,1,"1.pdb");

    // ---------- Create energy function
    let lj = PairwiseNonbondedEvaluator::new(1, CUTOFF as f32,
                                             Box::new(LennardJonesHomogenic::new(EPSILON_BY_K, SIGMA, CUTOFF)) );
    let mut total = TotalEnergy::default();
    total.add_component(Box::new(lj), 1.0);
    println!("{}", total.energy(&coords));

    // ---------- Create a sampler and add a mover into it
    let mut sampler = IsothermalMC::new(args.temperature as f64);
    sampler.energy = Box::new(total);               // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);

    // ---------- Run the simulation!
    let start = Instant::now();
    for i in 0..args.outer {
        let f_succ = sampler.run(&mut coords, args.inner as i32);
        to_pdb(&coords, i as i16, "ar.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&coords) / coords.size() as f64, f_succ, start.elapsed());
    }
}
