use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;

use crate::trait_acceptance_criterion::AcceptanceCriterion;


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
    pub fn new(temperature: f64) -> MetropolisCriterion {
        MetropolisCriterion {
            temperature,
            rng: SmallRng::from_entropy(),
        }
    }
}

impl AcceptanceCriterion for MetropolisCriterion {
    fn check(&mut self, energy_before: f64, energy_after: f64) -> bool {
        if energy_after <= energy_before {
            return true;
        } else {
            let delta_e = energy_after - energy_before;
            if self.rng.gen_range(0.0..1.0) < (-delta_e / self.temperature).exp() {
                return true;
            }
        }

        return false;
    }
}