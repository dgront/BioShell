
use std::ops::Range;
use std::default::Default;
use rand::{Rng, SeedableRng};
use rand::rngs::{SmallRng};

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

/// Metropolis acceptance criterion.
///
/// AcceptanceCriterion will return `true` or `false` when a Markov chain Monte Carlo move
/// from energy `energy_before` to `energy_after` should be accepted or not, respectively
pub trait AcceptanceCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool;
}

/// Classical Metropolis acceptance criterion.
///
/// A change of a system is accepted with probability `$P$`:
/// ```math
/// P(E_b \to E_a) = \begin{cases}\begin{align*}
///     1 &  \quad \text{when} \quad E_a \le E_b \\
///     e^{\Delta E / T} &  \quad \text{otherwise}
/// \end{align*}\end{cases}
/// ```
/// where `$\Delta E = E_a - E_b$` and the provided temperature `$T$` must be already expressed in
/// in the units of the Boltzmann constant `$k_B$`
#[derive(Clone)]
pub struct MetropolisCriterion {
    pub temperature: f64,
    rng: SmallRng,
}

impl MetropolisCriterion {
    /// Creates a new acceptance criterion that will result in Boltzmann distribution for the given temperature
    pub fn new(temperature: f64) -> MetropolisCriterion { MetropolisCriterion{temperature, rng:SmallRng::from_entropy() } }
}

impl AcceptanceCriterion for MetropolisCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool {

        if energy_after <= energy_before { return true; }
        else {
            let delta_e = energy_after - energy_before;
            if self.rng.gen_range(0.0..1.0) < (-delta_e / self.temperature).exp() { return true; }
        }

        return false;
    }
}

/// Mover changes a given conformation of a system
///
pub trait Mover<S: System, E: Energy<S>, T: AcceptanceCriterion> {

    /// Introduce a change into the given system
    ///
    /// The change may be accepted according to the acceptance criterion
    fn perturb(&mut self, system: &mut S, energy: &E, acc: &mut T) -> Option<Range<usize>>;

    /// The change may be accepted according to the acceptance criterion
    fn acceptance_statistics(&self) -> AcceptanceStatistics;

    /// Maximum range of perturbation allowed for that mover
    fn max_range(&self) -> f64;

    /// Sets the new maximum range of perturbation.
    fn set_max_range(&mut self, new_val: f64);
}

/// Defines a basic interface for a Monte Carlo sampling scheme.
pub trait Sampler<S: System, E: Energy<S>, T: AcceptanceCriterion> {

    /// Make `n` Monte Carlo sweeps
    ///
    /// During each sweep `coords.size()` random perturbations are attempted for each mover.
    /// Statistics for success rate are recorded in  [`AcceptanceStatistics`](AcceptanceStatistics)
    /// struct that is stored internally by each mover
    fn make_sweeps(&mut self, n:usize, coords: &mut S, energy: &E);

    /// Provides statistics of the success rate for this mover
    fn acceptance_criterion(&mut self) -> &mut T;

    /// Add a mover to this set
    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>);

    /// Immutable access to a mover
    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E, T>>;

    /// Mutable access to a mover
    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>>;

    /// Count movers contained in this set
    fn count_movers(&self) -> usize;

    fn acceptance_statistics(&self) -> Vec<AcceptanceStatistics> {

        let mut out:Vec<AcceptanceStatistics> = vec![];
        for i in 0..self.count_movers() {
            let stats = self.get_mover(i).acceptance_statistics().clone();
            out.push(stats);
        }
        return out;
    }
}


pub struct MCProtocol<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    acceptance_crit: T,
    movers: Vec<Box<dyn Mover<S, E, T>>>
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> MCProtocol<S, E, T> {
    pub fn new(acceptance_crit: T) -> MCProtocol<S, E, T> {
        MCProtocol {
            acceptance_crit, movers: vec![]
        }
    }
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> Sampler<S, E, T>  for MCProtocol<S, E, T> {
    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E) {
        for _ in 0..n {
            for i_mover in 0..self.movers.len() {
                let mover = &mut self.movers[i_mover];
                for _ in 0..coords.size() {
                    let _outcome = mover.perturb(coords, energy, &mut self.acceptance_crit);
                }
            }
        }
    }

    fn acceptance_criterion(&mut self) -> &mut T { &mut self.acceptance_crit }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>){
        self.movers.push(perturb_fn);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E, T>> { &self.movers[which_one] }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>> { &mut self.movers[which_one] }

    fn count_movers(&self) -> usize { self.movers.len() }
}


pub struct AdaptiveMCProtocol<S: System, E: Energy<S>, T: AcceptanceCriterion> {
    pub target_rate: f64,
    pub factor: f64,
    sampler: Box<dyn Sampler<S, E, T>>,
    allowed_ranges: Vec<Range<f64>>
}

impl<S: System, E: Energy<S>, T: AcceptanceCriterion> AdaptiveMCProtocol<S, E, T> {
    pub fn new(mut sampler: Box<dyn Sampler<S, E, T>>) -> AdaptiveMCProtocol<S, E, T> {
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

            let mut mover = self.sampler.get_mover_mut(i);
            let mut range = mover.max_range();
            if rate < self.target_rate - 0.05 { range = range * self.factor; }
            if rate > self.target_rate + 0.05 { range = range / self.factor; }
            if self.allowed_ranges[i].end.lt(&range) { range = self.allowed_ranges[i].end }
            if self.allowed_ranges[i].start.gt(&range) { range = self.allowed_ranges[i].start }
            mover.set_max_range(range);
        }
    }

    fn acceptance_criterion(&mut self) -> &mut T { self.sampler.acceptance_criterion() }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E, T>>) {
        let r = perturb_fn.max_range();
        self.sampler.add_mover(perturb_fn);
        self.allowed_ranges.push(r * 0.5..r * 4.0);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E, T>> { self.sampler.get_mover(which_one) }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E, T>> { self.sampler.get_mover_mut(which_one) }

    fn count_movers(&self) -> usize { self.sampler.count_movers() }
}
