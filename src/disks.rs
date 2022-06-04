use std::env;
use std::time::Instant;
use std::ops::Range;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb};
use bioshell_sim::generators::cubic_grid_atoms;
use bioshell_ff::nonbonded::SimpleContact;
use bioshell_sim::sampling::movers::{single_atom_move};
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

fn box_width(sphere_radius: f32, n_spheres: usize, density: f32) -> f32 {

    let v: f32 = 4.0/3.0 * std::f32::consts::PI * sphere_radius.powi(3);
    (n_spheres as f32 * v/density).powf(1.0/3.0)
}

pub fn main() {
    const R: f32 = 3.0;
    const W: f32 = 1.0;

    let n_atoms: usize = (6 as usize).pow(3);
    let density: f32 = 0.4;
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(R,n_atoms, density));
    cubic_grid_atoms(&mut coords);
    to_pdb(&coords,1,"1.pdb");

    // ---------- Create energy function
    let contacts = SimpleContact::new(R,R,R+W,1000.0,-1.0, 0);
    let mut total = TotalEnergy::default();
    total.add_component(Box::new(contacts), 1.0);
    println!("{}", total.energy(&coords));

    // ---------- Create a sampler and add a mover into it
    let mut sampler = IsothermalMC::new(0.7);
    sampler.energy = Box::new(total);               // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);

    // ---------- Run the simulation!
    let start = Instant::now();
    for i in 0..300 {
        let f_succ = sampler.run(&mut coords, 100);
        to_pdb(&coords, i, "disks.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&coords), f_succ, start.elapsed());
    }
}
