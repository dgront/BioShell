use std::ops::Range;

use bioshell_sim::{Energy, System};
use crate::acceptance_statistics::AcceptanceStatistics;
use crate::trait_mover::Mover;
use crate::trait_sampler::Sampler;

/// Adds to a given sampler the ability to adapt movers' range on the fly
///
/// This protocol monitors success rates for each mover contained in a given Monte Carlo [`Sampler`](Sampler)
/// and adjusts their `max_move_range` property to keep the success rate close to the desired value.
pub struct AdaptiveMCProtocol<S: System, E: Energy<S>> {
    pub target_rate: f64,
    pub factor: f64,
    sampler: Box<dyn Sampler<S, E>>,
    allowed_ranges: Vec<Range<f64>>,
}

impl<S: System, E: Energy<S>> AdaptiveMCProtocol<S, E> {
    pub fn new(sampler: Box<dyn Sampler<S, E>>) -> AdaptiveMCProtocol<S, E> {
        let mut allowed_ranges: Vec<Range<f64>> = vec![];
        for i in 0..sampler.count_movers() {
            let r = sampler.get_mover(i).max_range();
            allowed_ranges.push(r * 0.5..r * 4.0);
        }
        let out = AdaptiveMCProtocol {
            target_rate: 0.4,
            factor: 0.95,
            sampler,
            allowed_ranges,
        };
        return out;
    }
}

impl<S: System, E: Energy<S>> Sampler<S, E> for AdaptiveMCProtocol<S, E> {
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
            if rate < self.target_rate - 0.05 {
                range = range * self.factor;
            }
            if rate > self.target_rate + 0.05 {
                range = range / self.factor;
            }
            if self.allowed_ranges[i].end.lt(&range) {
                range = self.allowed_ranges[i].end
            }
            if self.allowed_ranges[i].start.gt(&range) {
                range = self.allowed_ranges[i].start
            }
            mover.set_max_range(range);
        }
    }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>) {
        let r = perturb_fn.max_range();
        self.sampler.add_mover(perturb_fn);
        self.allowed_ranges.push(r * 0.5..r * 4.0);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> {
        self.sampler.get_mover(which_one)
    }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> {
        self.sampler.get_mover_mut(which_one)
    }

    fn count_movers(&self) -> usize {
        self.sampler.count_movers()
    }
}