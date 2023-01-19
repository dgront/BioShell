
use std::ops::Range;
use std::default::Default;
use rand::Rng;
use rand::rngs::ThreadRng;

use bioshell_sim::{System, Energy};

#[derive(Clone, Debug)]
pub struct AcceptanceStatistics {
    pub n_succ:i32,
    pub n_failed:i32,
}

impl AcceptanceStatistics {
    pub fn success_rate(&self) -> f64 {
        let sum = self.n_succ + self.n_failed;
        if sum == 0 { return 0.0; }
        return self.n_succ as f64 / (sum as f64);
    }

    pub fn recent_success_rate(&self, prev_stats: &AcceptanceStatistics) -> f64 {
        let succ = self.n_succ - prev_stats.n_succ;
        let fail = self.n_failed - prev_stats.n_failed;
        let sum = succ + fail;
        if sum == 0 { return 0.0; }
        return succ as f64 / (sum as f64);
    }
}

impl Default for AcceptanceStatistics {
    fn default() -> Self {
        AcceptanceStatistics{ n_succ: 0, n_failed: 0 }
    }
}

pub trait AcceptanceCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool;
}

pub struct MetropolisCriterion {
    pub temperature: f64,
    rng: ThreadRng,
}

impl MetropolisCriterion {
    pub fn new(temperature: f64) -> MetropolisCriterion { MetropolisCriterion{temperature, rng:rand::thread_rng() } }
}

impl AcceptanceCriterion for MetropolisCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool {

        let delta_e = energy_after - energy_before;
        if delta_e <= 0.0 || self.rng.gen_range(0.0..1.0) <= (-delta_e / self.temperature).exp() {
            return true
        }
        return false;
    }
}

pub trait Mover<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    fn perturb(&mut self, system: &mut S, energy: &E, acc: &mut T) -> Option<Range<usize>>;
    fn acceptance_statistics(&self) -> AcceptanceStatistics;
    fn max_range(&self) -> f64;
    fn set_max_range(&mut self, new_val: f64);
}

pub trait Sampler<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    fn make_sweeps(&mut self, n:usize, coords: &mut S, energy: &E);
}

pub trait MoversSet<S: System, E: Energy<S>, T: AcceptanceCriterion> {

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>);

    fn get_mover(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>>;

    fn count_movers(&self) -> usize;
}

pub trait MoversSetSampler<S: System, E: Energy<S>, T: AcceptanceCriterion> : MoversSet<S, E, T> + Sampler<S, E, T> {}

pub struct MCProtocol<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    pub acceptance_criterion: T,
    movers: Vec<Box<dyn Mover<S, E, T>>>
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> MCProtocol<S, E, T> {
    pub fn new(acc_crit: T) -> MCProtocol<S, E, T> {
        MCProtocol {
            acceptance_criterion: acc_crit,
            movers: vec![]
        }
    }
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> Sampler<S, E, T>  for MCProtocol<S, E, T> {
    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E) {
        for _ in 0..n {
            for i_mover in 0..self.movers.len() {
                let mover = &mut self.movers[i_mover];
                for _ in 0..coords.size() {
                    let _outcome = mover.perturb(coords, energy, &mut self.acceptance_criterion);
                }
            }
        }
    }

    /*
    fn make_sweeps(&mut self, n:usize, coords: &mut S, energy: &Box<dyn Energy<S>>) {
        let mut future_coords = coords.clone();
        for _ in 0..n {
            for i_mover in 0..self.movers.len() {
                let mover = &mut self.movers[i_mover];
                for _ in 0..coords.size() {
                    // ---------- Make a move on future system
                    let range: Range<usize> = mover.perturb(&mut future_coords);
                    // ---------- Evaluate energy difference
                    let (en_before, en_after) = energy.delta_energy_by_range(coords, &future_coords, &range);
                    // ---------- test the energy consistency
                    #[cfg(debug_assertions)] {
                        let delta_en = en_after - en_before;
                        if delta_en < 1000.0 {      // ---------- test only when it's unlikely the move will be rejected
                            let en_old = energy.energy(&coords);
                            let en_new = energy.energy(&future_coords);
                            let total_delta = en_new - en_old;
                            if f64::abs(total_delta - delta_en) > 0.01 {
                                let str = format!("Inconsistent energy! Total {en_old} -> {en_new} with delta = {total_delta}, local delta: {delta_en} after move {:?}", range);
                                panic!("{}", str);
                            }
                        }
                    }
                    // ---------- apply acceptance criterion, copy or undo the move
                    if self.acceptance_criterion.check(en_before, en_after) {
                        // --- update mover counts, copy future_pose on current_pose to make the move
                        for ipos in range.start..range.end + 1 {
                            coords.copy_from(ipos, &future_coords);
                        }
                        mover.add_success();
                    } else {
                        // --- update mover failures, copy current_pose on future_pose to clear the move
                        for ipos in range.start..range.end + 1 {
                            future_coords.copy_from(ipos, &coords);
                        }
                        mover.add_failure();
                    }
                }
            }
        }
    }*/
}


impl<S: System, E: Energy<S>, T: AcceptanceCriterion> MoversSet<S, E, T> for MCProtocol<S, E, T> {

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>){
        self.movers.push(perturb_fn);
    }

    fn get_mover(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>> { &mut self.movers[which_one] }

    fn count_movers(&self) -> usize { self.movers.len() }
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> MoversSetSampler<S, E, T> for MCProtocol<S, E, T> {}

pub struct AdaptiveMCProtocol<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    pub target_rate: f64,
    pub factor: f64,
    sampler: Box<dyn MoversSetSampler<S, E, T>>,
    allowed_ranges: Vec<Range<f64>>
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> AdaptiveMCProtocol<S, E, T> {
    pub fn new(mut sampler: Box<dyn MoversSetSampler<S, E, T>>) -> AdaptiveMCProtocol<S, E, T> {
        let mut allowed_ranges:Vec<Range<f64>> = vec![];
        for i in 0..sampler.count_movers() {
            let r = sampler.get_mover(i).max_range();
            allowed_ranges.push(r * 0.5..r * 4.0);
        }
        let out = AdaptiveMCProtocol { target_rate: 0.4, factor: 0.95, sampler, allowed_ranges};
        return out;
    }
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> Sampler<S, E, T>  for AdaptiveMCProtocol<S, E, T> {

    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E) {
        let mut stats_before: Vec<AcceptanceStatistics> = vec![];
        for i in 0..self.sampler.count_movers() {
            stats_before.push(self.sampler.get_mover(i).acceptance_statistics());
        }
        self.sampler.make_sweeps(n, coords, energy);
        for i in 0..self.sampler.count_movers() {
            let stats_after = self.sampler.get_mover(i).acceptance_statistics();
            let rate = stats_after.recent_success_rate(&stats_before[i]);

            let mover = self.sampler.get_mover(i);
            let mut range = mover.max_range();
            if rate < self.target_rate - 0.05 { range = range * self.factor; }
            if rate > self.target_rate + 0.05 { range = range / self.factor; }
            if self.allowed_ranges[i].end.lt(&range) { range = self.allowed_ranges[i].end }
            if self.allowed_ranges[i].start.gt(&range) { range = self.allowed_ranges[i].start }
            mover.set_max_range(range);
        }
    }
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> MoversSet<S, E, T>  for AdaptiveMCProtocol<S, E, T> {
    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>) { self.sampler.add_mover(perturb_fn); }

    fn get_mover(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>> { self.sampler.get_mover(which_one) }

    fn count_movers(&self) -> usize { self.sampler.count_movers() }
}