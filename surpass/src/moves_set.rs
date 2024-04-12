use rand::Rng;
use rand::rngs::SmallRng;
use crate::{MoveProposal, SurpassAlphaSystem, SurpassEnergy, TotalEnergy};

struct MetropolisCriterion { temperature: f64}

impl MetropolisCriterion {
    pub fn new(temperature: f64) -> Self { MetropolisCriterion{temperature} }

    pub fn is_move_accepted(&self, energy_before: f64, energy_after: f64, rndgen: &mut SmallRng) -> bool {
        if energy_after < energy_before { return true }
        if rndgen.gen_range(0.0..1.0) < ((energy_before - energy_after) / self.temperature).exp() { return true; }
        return false;
    }
}

pub trait Mover {
    fn n_moved(&self) -> usize;
    fn propose(&self, system: &SurpassAlphaSystem, rndgen: &mut SmallRng, proposal: &mut MoveProposal);
}


pub struct MovesSet {
    movers: Vec<Box<dyn Mover>>,
    moves_in_cycle: Vec<usize>
}

impl MovesSet {
    pub fn new() -> Self { MovesSet{ movers: vec![], moves_in_cycle: vec![] } }

    pub fn add_mover(&mut self, energy_term: Box<dyn Mover>, moves_in_cycle: usize) {
        self.movers.push(energy_term);
        self.moves_in_cycle.push(moves_in_cycle);
    }

    pub fn mc_cycle(&self, system: &mut SurpassAlphaSystem, total_energy: &TotalEnergy, temperature: f64, rndgen: &mut SmallRng) {

        let criterion = MetropolisCriterion::new(temperature);

        for (mover, n_moves) in self.movers.iter().zip(&self.moves_in_cycle) {
            let mut proposal = MoveProposal::new(mover.n_moved());
            for _i in 0..*n_moves {
                mover.propose(system, rndgen, &mut proposal);
                let delta_en = total_energy.evaluate_delta(&system, &proposal);
                if criterion.is_move_accepted(0.0, delta_en, rndgen) {
                    proposal.apply(system);

                }
            }
        }
    }
}