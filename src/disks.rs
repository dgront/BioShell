use std::time::Instant;
use std::ops::Range;
use rand::Rng;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb, System};
use bioshell_ff::nonbonded::{NbList, PairwiseNonbondedEvaluator, SimpleContact, ArgonRules};
use bioshell_sim::generators::square_grid_atoms;
use bioshell_sim::sampling::protocols::{Ensemle, IsothermalMC, Sampler};

fn box_width(disc_radius: f64, n_discs: usize, density: f64) -> f64 {

    let v: f64 = std::f64::consts::PI * disc_radius.powi(2);
    (n_discs as f64 * v/density).powf(0.5)
}

/// Moves a random disc
fn single_dics_move(future: &mut System, max_step:f64) -> Range<usize> {
    let mut rng = rand::thread_rng();
    let i_moved = rng.gen_range(0..future.size());
    future.add(i_moved,rng.gen_range(-max_step..max_step),
               rng.gen_range(-max_step..max_step),0.0);

    i_moved..i_moved
}

pub fn main() {
    const R: f64 = 3.0;
    const W: f64 = 1.0;

    const N_SMALL_CYCLES: i32 = 100;
    const N_LARGE_CYCLES: i16 = 1000;

    let n_atoms: usize = 50 * 50;
    let density: f64 = 0.4;

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(R,n_atoms, density));
    square_grid_atoms(&mut coords);
    to_pdb(&coords,1,"disks.pdb");

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(R+W as f64,6.0,Box::new(ArgonRules{}));

    // ---------- Create the system
    let mut system: System = System::new(coords, nbl);

    // ---------- Create energy function
    let contacts = PairwiseNonbondedEvaluator::new(R+W,
        Box::new(SimpleContact::new(R,R,R+W,1000.0,-1.0)) );
    let mut total = TotalEnergy::default();
    total.add_component(Box::new(contacts), 1.0);
    println!("{} {} {}  {:.2?}", 0, total.energy(&system)/n_atoms as f64, 0.0, 0.0);

    // ---------- Create a sampler and add a mover into it
    let mut sampler = IsothermalMC::new(0.6, Ensemle::NVT, 1.0);
    sampler.energy = Box::new(total);               // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut System,f64) -> Range<usize>> = Box::new(single_dics_move);
    sampler.add_mover(m,3.0);

    // ---------- Run the simulation!
    let start = Instant::now();
    for i in 0..N_LARGE_CYCLES {
        let f_succ = sampler.run(&mut system, N_SMALL_CYCLES);
        to_pdb(&system.coordinates(), i, "disks.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&system)/system.size() as f64, f_succ, start.elapsed());
    }
}
