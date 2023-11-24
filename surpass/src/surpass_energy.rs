use crate::{MoveProposal, SurpassAlphaSystem};

/// Defines the energy function for the SURPASS-alpha model
pub trait SurpassEnergy {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64;
    fn evaluate_delta<const N_RES: usize, const N_ATOMS: usize>(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal<N_RES, N_ATOMS>) -> f64;
}

/// Defines the non-bonded energy kernel that evaluates energy between two atoms
pub trait NonBondedEnergyKernel {
    fn energy_for_distance_squared(&self, i2: f64) -> f64;

    fn distance_cutoff(&self) -> f64;
}