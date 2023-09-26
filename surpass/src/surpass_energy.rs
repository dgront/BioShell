use crate::{MoveProposal, SurpassAlphaSystem};

/// Defines the energy function for the SURPASS-alpha model
pub trait SurpassEnergy {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64;
    fn evaluate_delta<const N: usize>(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal<N>) -> f64;
}
