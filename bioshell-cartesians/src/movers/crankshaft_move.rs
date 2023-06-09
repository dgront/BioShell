use rand::Rng;
use std::ops::Range;
use crate::CartesianSystem;
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Energy, System};
use bioshell_numerical::{Rototranslation};

/// A mover that applies a crankshaft move to a molecule.
pub struct CrankshaftMove {
    frag_size: usize,
    max_angle: f64,
    succ_rate: AcceptanceStatistics,
}

impl CrankshaftMove {
    /// Create a new mover that applies a crankshaft move.
    pub fn new(max_angle: f64) -> CrankshaftMove {
        CrankshaftMove {
            max_angle,
            frag_size:5,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for CrankshaftMove {
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
        let mut rng = rand::thread_rng();
        let system_length = system.get_size();
        let i_start = rng.gen_range(0..system_length - self.frag_size - 1);
        let i_end = i_start + self.frag_size + 1;
        let angle = rng.gen_range(-self.max_angle..self.max_angle);
        let start = system.get_coordinates()[i_start].clone();
        let end = system.get_coordinates()[i_end].clone();

        // Take the periodic image of the end vector closest to the start vector
        let periodic_end = system.get_periodic_image(&start, &end);

        let roto_tran = Rototranslation::around_axis(&start, &periodic_end, angle);

        let energy_before = energy.energy(system);

        // Apply forward rotation
        for i in i_start + 1..i_end {
            let mut temp_coord = system.get_coordinates()[i].clone();
            roto_tran.apply_mut(&mut temp_coord);
            system.set_vec3(i, temp_coord);
        }

        let energy_after = energy.energy(system);

        if acc.check(energy_before, energy_after) {
            // Move succeeds
            self.succ_rate.n_succ += 1;
            Some(i_start..i_end + 1)
        } else {
            // Move fails
            self.succ_rate.n_failed += 1;
            // Fall back - apply inverse rotation
            for i in i_start + 1..i_end {
                let mut temp_coord = system.get_coordinates()[i].clone();
                roto_tran.apply_inverse_mut(&mut temp_coord);
                system.set_vec3(i, temp_coord);
            }
            None
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_angle
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_angle = new_val;
    }
}

