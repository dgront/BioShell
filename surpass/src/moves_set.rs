use log::debug;
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

    fn move_accepted(&mut self);

    fn move_cancelled(&mut self);

    fn adjust_move_range(&mut self);

    fn move_range(&mut self) -> f64;

    fn success_rate(&mut self) -> f64;

    fn reset_counters(&mut self) -> f64;
}

pub(crate) struct AdaptiveMoveRange {
    max_range: f64,
    max_range_allowed: f64,
    expected_acceptance: f64,
    n_accepted: usize,
    n_canceled: usize
}

impl AdaptiveMoveRange {
    pub fn new(max_range: f64, max_range_allowed: f64, expected_acceptance: f64) -> AdaptiveMoveRange {
        AdaptiveMoveRange {
            max_range, max_range_allowed,
            expected_acceptance,
            n_accepted: 0, n_canceled: 0,
        }
    }

    /// Record accepted move
    ///
    /// Statistics of how many moves were accepted and cancelled are used to adjust
    /// the maximum range of this mover.
    pub fn move_accepted(&mut self) { self.n_accepted += 1 }

    /// Record rejected move
    ///
    /// Statistics of how many moves were accepted and cancelled are used to adjust
    /// the maximum range of this mover.
    pub fn move_cancelled(&mut self) { self.n_canceled += 1 }

    pub fn reset_counters(&mut self) -> f64 {
        let out = self.success_rate();
        self.n_canceled = 0;
        self.n_accepted = 0;
        return out;
    }

    /// Adjust move range to reach the expected value of accepted moves
    pub fn adjust_move_range(&mut self) {
        let succ_rate = self.success_rate();

        if self.expected_acceptance - succ_rate > 0.05 {            // current success rate is too small, make moves smaller
            let tmp = self.max_range;
            self.max_range = self.max_range * 0.95;
            debug!("mover range shortened from {:5.2} to {:5.2} at acceptance rate {:6.3}",
                tmp, self.max_range, succ_rate);
        } else if succ_rate - self.expected_acceptance > 0.05 {     // current success rate is too high, make moves larger
            let tmp = self.max_range;
            self.max_range = self.max_range_allowed.min(self.max_range * 1.05);
            if self.max_range < self.max_range_allowed {
                debug!("mover range increased from {:5.2} to {:5.2} at acceptance rate {:6.3}",
                    tmp, self.max_range, succ_rate);
            } else {
                debug!("mover acceptance rate is {:5.3} but move range reached its allowed maximum {:6.3} range",
                    succ_rate, self.max_range_allowed);
            }
        }
    }

    pub fn move_range(&self) -> f64 { self.max_range }

    pub fn success_rate(&self) -> f64 { self.n_accepted as f64 / (self.n_accepted + self.n_canceled) as f64 }

}

pub struct MovesSet {
    movers: Vec<Box<dyn Mover>>,
    moves_in_cycle: Vec<usize>
}

impl MovesSet {
    pub fn new() -> Self { MovesSet{ movers: vec![], moves_in_cycle: vec![] } }

    pub fn add_mover(&mut self, mover: Box<dyn Mover>, moves_in_cycle: usize) {
        self.movers.push(mover);
        self.moves_in_cycle.push(moves_in_cycle);
    }

    pub fn mc_cycle(&mut self, system: &mut SurpassAlphaSystem, total_energy: &TotalEnergy,
        temperature: f64, n_cycles: usize, rndgen: &mut SmallRng) {
        let criterion = MetropolisCriterion::new(temperature);
        for _j in 0..n_cycles {
            for (mover, n_moves) in self.movers.iter_mut().zip(&self.moves_in_cycle) {
                let mut proposal = MoveProposal::new(mover.n_moved());
                for _i in 0..*n_moves {
                    mover.propose(system, rndgen, &mut proposal);
                    #[cfg(debug_assertions)]
                    check_bond_lengths(system, &proposal, 3.8);

                    let delta_en = total_energy.evaluate_delta(&system, &proposal);
                    if criterion.is_move_accepted(0.0, delta_en, rndgen) {
                        proposal.apply(system);
                        mover.move_accepted();
                    } else { mover.move_cancelled(); }
                }
            }
        }
        for mover in self.movers.iter_mut() {
            mover.adjust_move_range();
            mover.reset_counters();
        }
    }
}

#[allow(dead_code)]
fn check_bond_lengths(system: &mut SurpassAlphaSystem, mp: &MoveProposal, d: f64) {
    let mut backup: MoveProposal = MoveProposal::new(mp.n_moved);
    backup.first_moved_pos = mp.first_moved_pos;
    backup.backup(system);
    mp.apply(system);
    for i in 0..system.count_atoms()-1 {
        if system.chain(i) != system.chain(i+1) { continue }
        let dd = system.distance(i+1, i);
        if (dd-d).abs() > 0.01 {
            system.to_pdb_file("after.pdb", false);
            let av = system.ca_to_vec3(i+1);
            let prev_av = system.ca_to_vec3(i);
            backup.apply(system);
            let bv = system.ca_to_vec3(i+1);
            let prev_bv = system.ca_to_vec3(i);
            system.to_pdb_file("before.pdb", false);

            panic!("Broken bond between {} and {}, current length is: {}\nPos. before: {} {}\nPos. after: {} {}\n",
                   i, i+1, dd, &prev_bv, &bv, &prev_av, &av);
        }
    }
    backup.apply(system);
}

#[allow(dead_code)]
fn check_delta_en(en_before: f64, en_after: f64, delta: f64) {
    if (en_after-en_before-delta).abs() > 0.001 {
        panic!("Incorrect energy change: global {} vs delta {}\n", en_after-en_before, delta);
    }
}