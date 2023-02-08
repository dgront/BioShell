use log::{debug, info};

use std::ops::Range;
use std::default::Default;
use std::fmt;
use std::fmt::{Display};
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

/// Acceptance criterion for a Markov chain Monte Carlo
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
/// in the units of the Boltzmann constant `$k_B$`. This criterion, when used in atomistic
/// Monte Carlo simulation, results in NVT ensemble
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
pub trait Mover<S: System, E: Energy<S>> {

    /// Introduce a change into the given system
    ///
    /// The change may be accepted according to the acceptance criterion
    fn perturb(&mut self, system: &mut S, energy: &E, acc: &mut dyn AcceptanceCriterion) -> Option<Range<usize>>;

    /// The change may be accepted according to the acceptance criterion
    fn acceptance_statistics(&self) -> AcceptanceStatistics;

    /// Maximum range of perturbation allowed for that mover
    fn max_range(&self) -> f64;

    /// Sets the new maximum range of perturbation.
    fn set_max_range(&mut self, new_val: f64);
}

/// Defines a basic interface for a Monte Carlo sampling scheme.
pub trait Sampler<S: System, E: Energy<S>> {

    /// Make `n` Monte Carlo sweeps
    ///
    /// During each sweep `coords.size()` random perturbations are attempted for each mover added
    /// to this sampler.
    /// Statistics for success rate are recorded separately for each mover by
    /// [`AcceptanceStatistics`](AcceptanceStatistics) structs that are stored se by each mover.
    fn make_sweeps(&mut self, n:usize, coords: &mut S, energy: &E);

    /// Provides statistics of the success rate for this mover
    // fn acceptance_criterion(&mut self) -> &mut T;

    /// Add a mover to this set
    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>);

    /// Immutable access to a mover
    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>>;

    /// Mutable access to a mover
    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>>;

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


pub struct IsothermalMC<S: System, E: Energy<S>> {
    acceptance_crit: MetropolisCriterion,
    movers: Vec<Box<dyn Mover<S, E>>>
}

impl<S: System, E: Energy<S>> IsothermalMC<S, E> {
    pub fn new(temperature: f64) -> IsothermalMC<S, E> {
        IsothermalMC {
            acceptance_crit: MetropolisCriterion::new(temperature),
            movers: vec![]
        }
    }
}

impl<S: System, E: Energy<S>> Sampler<S, E>  for IsothermalMC<S, E> {
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

    // fn acceptance_criterion(&mut self) -> &mut T { &mut self.acceptance_crit }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>){
        self.movers.push(perturb_fn);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> { &self.movers[which_one] }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> { &mut self.movers[which_one] }

    fn count_movers(&self) -> usize { self.movers.len() }
}


pub struct AdaptiveMCProtocol<S: System, E: Energy<S>> {
    pub target_rate: f64,
    pub factor: f64,
    sampler: Box<dyn Sampler<S, E>>,
    allowed_ranges: Vec<Range<f64>>
}

impl<S: System, E: Energy<S>> AdaptiveMCProtocol<S, E> {
    pub fn new(sampler: Box<dyn Sampler<S, E>>) -> AdaptiveMCProtocol<S, E> {
        let mut allowed_ranges:Vec<Range<f64>> = vec![];
        for i in 0..sampler.count_movers() {
            let r = sampler.get_mover(i).max_range();
            allowed_ranges.push(r * 0.5..r * 4.0);
        }
        let out = AdaptiveMCProtocol { target_rate: 0.4, factor: 0.95, sampler, allowed_ranges};
        return out;
    }
}

impl<S: System, E: Energy<S>> Sampler<S, E>  for AdaptiveMCProtocol<S, E> {

    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E) {
        let mut stats_before: Vec<AcceptanceStatistics> = vec![];
        for i in 0..self.sampler.count_movers() {
            stats_before.push(self.sampler.get_mover(i).acceptance_statistics());
        }
        self.sampler.make_sweeps(n, coords, energy);
        for i in 0..self.sampler.count_movers() {
            let stats_after = self.sampler.get_mover(i).acceptance_statistics();
            let rate = stats_after.recent_success_rate(&stats_before[i]);

            let mover = self.sampler.get_mover_mut(i);
            let mut range = mover.max_range();
            if rate < self.target_rate - 0.05 { range = range * self.factor; }
            if rate > self.target_rate + 0.05 { range = range / self.factor; }
            if self.allowed_ranges[i].end.lt(&range) { range = self.allowed_ranges[i].end }
            if self.allowed_ranges[i].start.gt(&range) { range = self.allowed_ranges[i].start }
            mover.set_max_range(range);
        }
    }

    // fn acceptance_criterion(&mut self) -> &mut T { self.sampler.acceptance_criterion() }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>) {
        let r = perturb_fn.max_range();
        self.sampler.add_mover(perturb_fn);
        self.allowed_ranges.push(r * 0.5..r * 4.0);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> { self.sampler.get_mover(which_one) }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> { self.sampler.get_mover_mut(which_one) }

    fn count_movers(&self) -> usize { self.sampler.count_movers() }
}

pub trait StepwiseMover<S: System, E: Energy<S>> {
    fn start(&mut self, system: &mut S, energy: &E) -> f64;
    fn grow_by_one(&mut self, system: &mut S, energy: &E) -> f64;
}

pub trait StepwiseBuilder<S: System, E: Energy<S>> {

    fn build(&mut self, system: &mut S, energy: &E) -> f64;
}

#[allow(non_snake_case)]
pub struct PERM<S: System, E: Energy<S>> {
    /// controls the pruning rate
    pub c_low: f64,
    /// controls the enrichment rate
    pub c_hi: f64,
    /// the number of full chains built do far
    chains_cnt: i64,
    Z: Vec<f64>,                        // --- partition function for each chain length
    W_low: Vec<f64>,                    // --- lower bound for the W (pruning criteria)
    W_hi: Vec<f64>,                     // --- upper bound for the W (enrichment criteria)
    step: Box<dyn StepwiseMover<S, E>>, // --- stepwise mover used to build up chains
    chains: Vec<(S, f64)>               // --- stack for enriched chains
}

impl<S: System, E: Energy<S>> PERM<S, E> {

    #[allow(non_snake_case)]
    pub fn new(n:usize, c_low: f64, c_hi: f64, step: Box<dyn StepwiseMover<S, E>>) -> PERM<S, E> {
        let Z = vec![0.0; n];
        let W_low = vec![0.0; n];
        let large = 1000000.0;
        let mut W_hi = vec![large; n];
        for i in 1..W_hi.len() { W_hi[i] = W_hi[i - 1] * large; }
        let chains = Vec::new();
        return PERM{c_low, c_hi, chains_cnt: 0, Z, W_low, W_hi, chains, step};
    }

    /// Returns the capacity (maximum size) of system generated by this PERM generator
    pub fn capacity(&self) -> usize { self.W_low.len() }

    /// Returns the number of chains created so far
    pub fn count_chains(&self) -> i64 { self.chains_cnt }

    /// Updates the internal weights that control pruning and enrichment events
    pub fn update_weights(&mut self) {
        for i in 0..self.capacity() {
            self.W_low[i] = self.c_low * self.Z[i] / self.Z[0];
            self.W_hi[i] = self.c_hi * self.Z[i] / self.Z[0];
        }
    }

    pub fn chains_left(&self) -> bool { !self.chains.is_empty() }
}

impl<S: System, E: Energy<S>> Display for PERM<S, E> {
    /// Creates a `String` representation of this `PERM` sampler.
    /// The output shows the current state of the internal data: i, Z, W_lo, W_hi in respective columns
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::new();
        for i in 0..self.capacity() {
            out = format!("{}{:.6} {:.4} {:.2} {:.2} {:.2}\n", out,
                          self.chains_cnt, i, self.Z[i], self.W_low[i], self.W_hi[i]);
        }
        write!(f, "{}", out)
    }
}

impl<S: System, E: Energy<S>> StepwiseBuilder<S, E> for PERM<S, E> {

    /// Generate a new system
    ///
    /// The generated state (e.g. new coordinates) will be stored in the given `system` reference;
    /// statistical weight of this new chain will be returned.
    fn build(&mut self, system: &mut S, energy: &E) -> f64 {
        let mut current_pos: usize;     // --- grow the current chain from the current_pos to the maximum length (i.e. capacity)
        let mut w_tot: f64;             // --- total (cumulative) weight of the current chain
        // ---------- take the most recent stub from the stack, if there is one there waiting
        if self.chains.len() > 0 {
            let (sys, mut w) = self.chains.pop().unwrap();
            for i in 0..sys.size() {
                system.copy_from(i, &sys);
            }
            current_pos = sys.size();
            system.set_size(current_pos);
            w_tot = w;
        } else {
            current_pos = 1;
            self.step.start(system, energy);
            w_tot = 1.0;
            self.Z[0] += 1.0;
            debug!("starting a new system");
        }

        for i in current_pos..self.capacity() {
            // ---------- add next element of the system
            let w: f64 = self.step.grow_by_one(system, energy);
            w_tot *= w;
            if w_tot < self.W_low[i] {      // --- prune event
                let mut rng = rand::thread_rng();
                let if_prune = rng.gen_range(0.0..1.0);
                if if_prune < 0.5 { return 0.0;}
                else { w_tot *= 2.0; }
            }
            if w_tot > self.W_hi[i] {      // --- enrich event
                w_tot *= 0.5;
                self.chains.push((system.clone(), w_tot));
                info!("enrichment: pushing on stack a cloned system of size {}",system.size());
            }
            self.Z[i] += w_tot;
        }

        return w_tot;
    }
}






