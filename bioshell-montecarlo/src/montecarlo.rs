use log::{debug, info};

use std::ops::Range;
use std::default::Default;
use std::fmt;
use std::fmt::{Display};
use std::time::Instant;
use rand::{Rng, SeedableRng};
use rand::rngs::{SmallRng};

use bioshell_sim::{System, Energy, ResizableSystem, ObserversSet};

#[derive(Clone, Debug)]
/// Counts how many system perturbations were successful.
///
/// Each Monte Carlo [`Mover`](Mover) must contain an [AcceptanceStatistics]
/// and update its counters accordingly to the outcome of a [`Mover::perturb()`] call.
/// The total number of Monte Carlo moves attempted is `n_succ + n_failed`
pub struct AcceptanceStatistics {
    /// number of successful perturbations
    pub n_succ:i32,
    /// number of failures
    pub n_failed:i32,
}

impl AcceptanceStatistics {
    /// Computes the success rate for a given Monte Carlo Markov chain.
    ///
    /// Simply returns `n_succ / (n_succ + n_failed)`
    pub fn success_rate(&self) -> f64 {
        let sum = self.n_succ + self.n_failed;
        if sum == 0 { return 0.0; }
        return self.n_succ as f64 / (sum as f64);
    }

    /// Computes the success rate since the given point in simulation
    ///
    /// The success rate is computed based on *new* observations that were made
    /// after the given `prev_stats` were recorded
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

    /// Run a Monte Carlo simulation
    ///
    /// A simulation consists of a `n_outer` rounds of `n_inner` Monte Carlo sweeps. After each
    /// round of Monte Carlo sweeps observations are taken.
    fn run_simulation(&mut self, n_inner:usize, n_outer: usize, coords: &mut S, energy: &E, observers: &mut ObserversSet<S>) {
        // ---------- Run the simulation!
        let start = Instant::now();
        let mut recent_acceptance: Vec<AcceptanceStatistics> = vec![AcceptanceStatistics::default(); self.count_movers()];
        for i in 0..n_outer {
            self.make_sweeps(n_inner, coords, &energy);
            print!("{:6} {:9.3}  ", i, energy.energy(&coords) / coords.size() as f64);
            for i_mover in 0..self.count_movers() {
                let stats = self.get_mover(i_mover).acceptance_statistics();
                print!("{:5.3} ", stats.recent_success_rate(&recent_acceptance[i_mover]));
                recent_acceptance[i_mover] = stats;
            }
            println!(" {:.2?}", start.elapsed());
            observers.observe(coords);
        }
    }

    /// Add a mover to this set
    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>);

    /// Immutable access to a mover
    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>>;

    /// Mutable access to a mover
    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>>;

    /// Count movers contained in this set
    fn count_movers(&self) -> usize;
}

/// Isothermal Monte Carlo simulation.
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

    /// Returns temperature of this simulation
    pub fn temperature(&self) -> f64 { self.acceptance_crit.temperature }
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

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>){
        self.movers.push(perturb_fn);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> { &self.movers[which_one] }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> { &mut self.movers[which_one] }

    fn count_movers(&self) -> usize { self.movers.len() }
}


/// Adds to a given sampler the ability to adapt movers' range on the fly
///
/// This protocol monitors success rates for each mover contained in a given Monte Carlo [`Sampler`](Sampler)
/// and adjusts their `max_move_range` property to keep the success rate close to the desired value.
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

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>) {
        let r = perturb_fn.max_range();
        self.sampler.add_mover(perturb_fn);
        self.allowed_ranges.push(r * 0.5..r * 4.0);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> { self.sampler.get_mover(which_one) }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> { self.sampler.get_mover_mut(which_one) }

    fn count_movers(&self) -> usize { self.sampler.count_movers() }
}

pub trait StepwiseMover<S: ResizableSystem, E: Energy<S>> {
    fn start(&mut self, system: &mut S, energy: &E) -> f64;
    fn grow_by_one(&mut self, system: &mut S, energy: &E) -> f64;
}

pub trait StepwiseBuilder<S: ResizableSystem, E: Energy<S>> {

    fn build(&mut self, system: &mut S, energy: &E) -> f64;
}

#[allow(non_snake_case)]
pub struct PERM<S: ResizableSystem, E: Energy<S>> {
    /// controls the pruning rate
    pub c_low: f64,
    /// controls the enrichment rate
    pub c_hi: f64,
    /// [update_weights()]`update_weights()` method is called every `n` full chains generated
    pub update_after_n: usize,
    /// the number of full chains built do far
    chains_cnt: i64,
    Z: Vec<f64>,                        // --- partition function for each chain length
    W_low: Vec<f64>,                    // --- lower bound for the W (pruning criteria)
    W_hi: Vec<f64>,                     // --- upper bound for the W (enrichment criteria)
    Z_div_step: f64,                    // --- at every k-th step, the Z[k] value is divided by Z_div_step^k to make sure it fits info f64
    step: Box<dyn StepwiseMover<S, E>>, // --- stepwise mover used to build up chains
    chains: Vec<(S, f64)>               // --- stack for enriched chains
}

impl<S: ResizableSystem, E: Energy<S>> PERM<S, E> {

    #[allow(non_snake_case)]
    pub fn new(n:usize, c_low: f64, c_hi: f64, step: Box<dyn StepwiseMover<S, E>>) -> PERM<S, E> {
        // ---------- start partition function with zeros
        let Z = vec![0.0; n];
        // ---------- start pruning weights with zeros
        let W_low = vec![0.0; n];
        // ---------- start enriching weights with something extremely large
        let mut W_hi = vec![1000000.0; n];
        // ---------- Z_div[k] =  10^k
        let Z_div_step: f64 = 10.0;
        for i in 1..W_hi.len() {
            W_hi[i] = W_hi[i - 1] * 1000000.0;
        }
        let chains = Vec::new();
        return PERM{c_low, c_hi, update_after_n:1, chains_cnt: 0, Z, W_low, W_hi, Z_div_step, chains, step};
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

    /// Provides statistics for the `i`-th step of generated chains.
    ///
    /// Returns `$W_i^L$`, `$W_i^H$` and `$Z_i / Z_0$`, i.e. the two critical values for pruning
    /// and enriching events (`$W_i^L$` and `$W_i^H$`, respectively) as well as the current value
    /// of the (normalized) partition function.
    pub fn weights(&self, i: usize) -> (f64, f64, f64) {
        (self.W_low[i], self.W_hi[i], self.Z[i] / self.Z[0])
    }

    pub fn chains_left(&self) -> bool { !self.chains.is_empty() }

    /// Sets a value each `w[i]` is divided by.
    ///
    /// The [`PERM`](PERM) method accumulates the total weight for a growing system as a product of weights
    /// obtained at each step, i.e. by adding a new element to a chain. Even for a modest value of
    /// a per-element weight, their product for a long chain can exceed the capacity of the `f64` variable.
    /// Therefore, the weight for a chain of `k` elements is stored internally by [`PERM`](PERM)
    /// as:
    /// ```math
    /// w = \prod_1^{k} \frac{w_i}{c} = \frac{1}{c^k} \prod_1^{k} \frac{w_i}
    /// ```
    /// By default, `c=10`; however it's reasonable to set the `c` constant to the expected value of
    /// `$w_i$`
    pub fn set_w_scale(&mut self, w_div: f64) { self.Z_div_step = w_div; }
}

impl<S: ResizableSystem, E: Energy<S>> Display for PERM<S, E> {
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

impl<S: ResizableSystem, E: Energy<S>> StepwiseBuilder<S, E> for PERM<S, E> {

    /// Generate a new system
    ///
    /// The generated state (e.g. new coordinates) will be stored in the given `system` reference;
    /// statistical weight of this new chain will be returned.
    ///
    /// If the current chain has been pruned, this method returns 0.0 for its weight.
    fn build(&mut self, system: &mut S, energy: &E) -> f64 {
        let current_pos: usize;         // --- grow the current chain from the current_pos to the maximum length (i.e. capacity)
        let mut w_tot: f64;             // --- total (cumulative) weight of the current chain
        // ---------- take the most recent stub from the stack, if there is one there waiting
        if self.chains.len() > 0 {
            let (sys, w) = self.chains.pop().unwrap();
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
            // ---------- add the next element of the system
            let w: f64 = self.step.grow_by_one(system, energy);
            w_tot *= w / self.Z_div_step;                                                         // --- adjust the weight so it fits info f64
            if w_tot < self.W_low[i] && self.chains_cnt > (self.update_after_n * 2) as i64 {      // --- prune event
                let mut rng = rand::thread_rng();
                let if_prune = rng.gen_range(0.0..1.0);
                if if_prune < 0.5 { return 0.0;}
                else { w_tot *= 2.0; }
            }
            if w_tot > self.W_hi[i]  && self.chains_cnt > (self.update_after_n * 2) as i64 {      // --- enrich event
                w_tot *= 0.5;
                self.chains.push((system.clone(), w_tot));
                debug!("enrichment: pushing on stack a cloned system of size {}",system.size());
            }
            self.Z[i] += w_tot;
        }
        if self.chains_cnt % self.update_after_n as i64 == (self.update_after_n - 1) as i64 {
            self.update_weights();
        }
        self.chains_cnt += 1;
        return w_tot;
    }
}






