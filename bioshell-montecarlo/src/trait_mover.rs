use std::ops::Range;

use bioshell_sim::{Energy, System};
use crate::trait_acceptance_criterion::AcceptanceCriterion;
use crate::acceptance_statistics::AcceptanceStatistics;


/// Mover changes a given conformation of a system
///
pub trait Mover<S: System, E: Energy<S>> {
    /// Introduce a change into the given system
    ///
    /// The change may be accepted according to the acceptance criterion
    fn perturb(
        &mut self,
        system: &mut S,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>>;

    /// The change may be accepted according to the acceptance criterion
    fn acceptance_statistics(&self) -> AcceptanceStatistics;

    /// Maximum range of perturbation allowed for that mover
    fn max_range(&self) -> f64;

    /// Sets the new maximum range of perturbation.
    fn set_max_range(&mut self, new_val: f64);
}