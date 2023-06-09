/// Acceptance criterion for a Markov chain Monte Carlo
///
/// AcceptanceCriterion will return `true` or `false` when a Markov chain Monte Carlo move
/// from energy `energy_before` to `energy_after` should be accepted or not, respectively
pub trait AcceptanceCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool;
}