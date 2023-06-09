use std::default::Default;
use std::time::Instant;

use bioshell_sim::{Energy, ObserversSet, System};
use crate::acceptance_statistics::AcceptanceStatistics;
use crate::trait_mover::Mover;

/// Defines a basic interface for a Monte Carlo sampling scheme.
pub trait Sampler<S: System, E: Energy<S>> {
    /// Make `n` Monte Carlo sweeps
    ///
    /// During each sweep `coords.size()` random perturbations are attempted for each mover added
    /// to this sampler.
    /// Statistics for success rate are recorded separately for each mover by
    /// [`AcceptanceStatistics`](AcceptanceStatistics) structs that are stored se by each mover.
    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E);

    /// Run a Monte Carlo simulation
    ///
    /// A simulation consists of a `n_outer` rounds of `n_inner` Monte Carlo sweeps. After each
    /// round of Monte Carlo sweeps observations are taken.
    fn run_simulation(
        &mut self,
        n_inner: usize,
        n_outer: usize,
        coords: &mut S,
        energy: &E,
        observers: &mut ObserversSet<S>,
    ) {
        // ---------- Run the simulation!
        let start = Instant::now();
        let mut recent_acceptance: Vec<AcceptanceStatistics> =
            vec![AcceptanceStatistics::default(); self.count_movers()];
        for i in 0..n_outer {
            self.make_sweeps(n_inner, coords, &energy);
            print!(
                "{:6} {:9.3}  ",
                i,
                energy.energy(&coords) / coords.get_size() as f64
            );
            for i_mover in 0..self.count_movers() {
                let stats = self.get_mover(i_mover).acceptance_statistics();
                print!(
                    "{:5.3} ",
                    stats.recent_success_rate(&recent_acceptance[i_mover])
                );
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