use bioshell_sim::{Energy, System};
use crate::metropolis_criterion::MetropolisCriterion;
use crate::trait_mover::Mover;
use crate::trait_sampler::Sampler;

/// Isothermal Monte Carlo simulation.
pub struct IsothermalMC<S: System, E: Energy<S>> {
    acceptance_crit: MetropolisCriterion,
    movers: Vec<Box<dyn Mover<S, E>>>,
}

impl<S: System, E: Energy<S>> IsothermalMC<S, E> {
    pub fn new(temperature: f64) -> IsothermalMC<S, E> {
        IsothermalMC {
            acceptance_crit: MetropolisCriterion::new(temperature),
            movers: vec![],
        }
    }

    /// Returns temperature of this simulation
    pub fn temperature(&self) -> f64 {
        self.acceptance_crit.temperature
    }
}

impl<S: System, E: Energy<S>> Sampler<S, E> for IsothermalMC<S, E> {
    fn make_sweeps(&mut self, n: usize, coords: &mut S, energy: &E) {
        for _ in 0..n {
            for i_mover in 0..self.movers.len() {
                let mover = &mut self.movers[i_mover];
                for _ in 0..coords.get_size() {
                    let _outcome = mover.perturb(coords, energy, &mut self.acceptance_crit);
                }
            }
        }
    }

    fn add_mover(&mut self, perturb_fn: Box<dyn Mover<S, E>>) {
        self.movers.push(perturb_fn);
    }

    fn get_mover(&self, which_one: usize) -> &Box<dyn Mover<S, E>> {
        &self.movers[which_one]
    }

    fn get_mover_mut(&mut self, which_one: usize) -> &mut Box<dyn Mover<S, E>> {
        &mut self.movers[which_one]
    }

    fn count_movers(&self) -> usize {
        self.movers.len()
    }
}