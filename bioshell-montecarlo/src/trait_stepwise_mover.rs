//use crate::trait_stepwise_builder::StepwiseBuilder;
use bioshell_sim::{Energy, ResizableSystem};

pub trait StepwiseMover<S: ResizableSystem, E: Energy<S>> {
    fn start(&mut self, system: &mut S, energy: &E) -> f64;
    fn grow_by_one(&mut self, system: &mut S, energy: &E) -> f64;
}