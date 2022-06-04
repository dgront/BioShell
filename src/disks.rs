use std::time::Instant;
use std::ops::Range;
use rand::Rng;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb};
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_sim::generators::square_grid_atoms;
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

fn box_width(disc_radius: f32, n_discs: usize, density: f32) -> f32 {

    let v: f32 = std::f32::consts::PI * disc_radius.powi(2);
    (n_discs as f32 * v/density).powf(0.5)
}

/// Moves a random disc
fn single_dics_move(future: &mut Coordinates, max_step:f32) -> Range<usize> {
    let mut rng = rand::thread_rng();
    let i_moved = rng.gen_range(0..future.size());
    future.add(i_moved,rng.gen_range(-max_step..max_step),
               rng.gen_range(-max_step..max_step),0.0);

    i_moved..i_moved
}

pub fn main() {
    const R: f32 = 3.0;
    const W: f32 = 1.0;

    let n_atoms: usize = 50 * 50;
    let density: f32 = 0.4;
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(R,n_atoms, density));
    square_grid_atoms(&mut coords);
    to_pdb(&coords,1,"1.pdb");

    // ---------- Create energy function
    let contacts = PairwiseNonbondedEvaluator::new(1,R+W,
        Box::new(SimpleContact::new(R,R,R+W,1000.0,-1.0)) );
    let mut total = TotalEnergy::default();
    total.add_component(Box::new(contacts), 1.0);
    println!("{}", total.energy(&coords));

    // ---------- Create a sampler and add a mover into it
    let mut sampler = IsothermalMC::new(0.6);
    sampler.energy = Box::new(total);               // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(single_dics_move);
    sampler.add_mover(m,3.0);

    // ---------- Run the simulation!
    let start = Instant::now();
    for i in 0..300 {
        let f_succ = sampler.run(&mut coords, 100);
        to_pdb(&coords, i, "disks.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&coords)/coords.size() as f64, f_succ, start.elapsed());
    }
}
