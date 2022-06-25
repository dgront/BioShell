use std::env;
use std::time::Instant;
use std::ops::Range;

use bioshell_numerical::Vec3;
use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb, System};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::{SimpleContact, PairwiseNonbondedEvaluator, NbList, NbListRules};
use bioshell_sim::generators::random_chain;
use bioshell_sim::sampling::movers::{single_atom_move, perturb_chain_fragment};
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

#[derive(Clone)]
pub struct PolymerRules;

impl NbListRules for PolymerRules {

    fn if_atom_excluded(&self, _coordinates: &Coordinates, _i_atom: usize) -> bool { false }

    fn if_pair_excluded(&self, _coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool {
        if i_atom >  j_atom { i_atom - j_atom > 2} else { j_atom - i_atom > 2 }
    }

    fn box_clone(&self) -> Box<dyn NbListRules> {
        Box::new((*self).clone())
    }
}

pub fn main() {
    const E_REP: f32 = 4.0;
    const E_FROM: f32 = 4.5;
    const E_TO: f32 = 6.0;
    const E_VAL: f32 = -1.0;
    const REP_VAL: f32 = 1000.0;
    const N: usize = 30;

    let args: Vec<String> = env::args().collect();

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(N);
    coords.set_box_len(1000000.0);
    let start: Vec3 = Vec3::new(100.0,100.0,100.0);
    random_chain(3.8, E_REP as f64, &start, &mut coords);
    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(E_TO,3.0,Box::new(PolymerRules{}));
    // ---------- Create the system
    let mut system: System = System::new(coords, nbl);


    // chain_to_pdb(&coords,"1.pdb");

    // ---------- Contact energy
    let contacts = PairwiseNonbondedEvaluator::new(1, 6.0 as f32,
            Box::new(SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL)));

    let harmonic = SimpleHarmonic::new(3.8,2.0);


    let mut total = TotalEnergy::default();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    println!("{}", total.energy(&system));

    let mut sampler = IsothermalMC::new(2.0);
    sampler.energy = Box::new(total);    // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut System,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);
    // let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(perturb_chain_fragment);
    // sampler.add_mover(m,5.0);

    let start = Instant::now();
    for i in 0..1000 {
        let f_succ = sampler.run(&mut system, 100);
        to_pdb(&system.coordinates(), i, "tra.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&system), f_succ, start.elapsed());
    }
}
