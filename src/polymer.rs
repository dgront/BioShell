use std::time::Instant;
use std::ops::Range;

use bioshell_numerical::Vec3;
use bioshell_ff::{Coordinates, Energy, TotalEnergy, to_pdb, System};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::{SimpleContact, PairwiseNonbondedEvaluator, NbList, NbListRules};
use bioshell_sim::generators::random_chain;
use bioshell_sim::sampling::movers::{single_atom_move};
use bioshell_sim::sampling::protocols::{Ensemle, IsothermalMC, Sampler};

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
    const E_REP: f32 = 4.25;
    const E_FROM: f32 = 4.5;
    const E_TO: f32 = 6.0;
    const E_VAL: f32 = -1.0;
    const REP_VAL: f32 = 1000.0;
    const N_BEADS: usize = 300;
    const N_SMALL_CYCLES: i32 = 1000;
    const N_LARGE_CYCLES: i16 = 10000;
    const T: f64 = 1.8;
    const L: f32 = 900.0;

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(N_BEADS);
    coords.set_box_len(L);
    let start: Vec3 = Vec3::new(L/2.0, L/2.0, L/2.0);
    random_chain(3.8, E_FROM as f64, &start, &mut coords);
    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(E_TO,4.0,Box::new(PolymerRules{}));
    // ---------- Create the system
    let mut system: System = System::new(coords, nbl);

    to_pdb(&system.coordinates(), 0, "tra.pdb");

    // ---------- Contact energy
    let contacts = PairwiseNonbondedEvaluator::new(E_TO as f32,
            Box::new(SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL)));
    // ---------- Harmonic energy (i.e. springs between beads)
    let harmonic = SimpleHarmonic::new(3.8,5.0);

    let mut total = TotalEnergy::default();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    println!("{}", total.energy(&system));

    let mut sampler = IsothermalMC::new(T, Ensemle::NVT, 1.0);
    sampler.energy = Box::new(total);    // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut System,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);
    // let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(perturb_chain_fragment);
    // sampler.add_mover(m,5.0);

    let start = Instant::now();
    for i in 0..N_LARGE_CYCLES {
        let f_succ = sampler.run(&mut system, N_SMALL_CYCLES);
        to_pdb(&system.coordinates(), i+1, "tra.pdb");
        println!("{} {} {}  {:.2?}", i, sampler.energy(&system), f_succ, start.elapsed());
    }
}
