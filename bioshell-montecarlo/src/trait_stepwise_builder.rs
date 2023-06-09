use bioshell_sim::{Energy, ResizableSystem};

pub trait StepwiseBuilder<S: ResizableSystem, E: Energy<S>> {
    fn build(&mut self, system: &mut S, energy: &E) -> f64;
}