use rand::Rng;
use std::ops::Range;
use crate::CartesianSystem;
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Energy, System};

/// A mover that moves a single, randomly selected atom by a small random vector.
pub struct SingleAtomMove {
    max_step: f64,
    succ_rate: AcceptanceStatistics,
}

impl SingleAtomMove {
    /// Create a new mover that shifts a single atom.
    pub fn new(max_range: f64) -> SingleAtomMove {
        SingleAtomMove {
            max_step: max_range,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for SingleAtomMove {
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
        let mut rng = rand::thread_rng();//obtain a random number generator
        let i_moved = rng.gen_range(0..system.get_size());//obtain a random number

        let old_x: f64 = system.get_coordinates()[i_moved].x;//obtain a random coordinate of x
        let old_y: f64 = system.get_coordinates()[i_moved].y;//obtain a random coordinate of y
        let old_z: f64 = system.get_coordinates()[i_moved].z;//obtain a random coordinate of z
        let old_energy: f64 = energy.energy_by_pos(system, i_moved);

        //add the randomly generated coordinate into the Cartesian System.
        system.add_xyz(
            i_moved,
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
        );

        let new_energy: f64 = energy.energy_by_pos(system, i_moved);

        if acc.check(old_energy, new_energy)
        {
            //success
            system.update_nbl(i_moved);//update non-bonded list of neighbors
            self.succ_rate.n_succ += 1;
            return Option::from(i_moved..i_moved);
        }
        else
        {
            //failure
            self.succ_rate.n_failed += 1;
            system.set_xyz(i_moved, old_x, old_y, old_z);//fall back
            return Option::None;
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_step
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_step = new_val;
    }
}