use std::env;
use std::time::Instant;
use std::ops::Range;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::SimpleContact;
use bioshell_sim::generators::random_chain;
use bioshell_sim::sampling::movers::{single_atom_move, perturb_chain_fragment};
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

pub fn main() {
    let args: Vec<String> = env::args().collect();

    let n_res: usize = 300;
    let mut coords = Coordinates::new(n_res);
    random_chain(3.8,4.5,&mut coords);
    // chain_to_pdb(&coords,"1.pdb");

    let harmonic = SimpleHarmonic::new(3.8,2.0);
    let contacts = SimpleContact::new(4.0,4.5,6.0,1000.0,-1.0, 1);

    let mut total = TotalEnergy::default();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    println!("{}", total.energy(&coords));

    let mut sampler = IsothermalMC::new(2.0);
    sampler.energy = Box::new(total);    // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);
    // let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(perturb_chain_fragment);
    // sampler.add_mover(m,5.0);

    let start = Instant::now();
    for i in 0..1000 {
        let f_succ = sampler.run(&mut coords, 100);
        to_pdb(&coords, i, "tra.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&coords), f_succ, start.elapsed());
    }
}
