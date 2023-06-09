//! Provides means to build a Cartesian system, either in a deterministic or a stochastic way
use bioshell_montecarlo::{StepwiseBuilder, StepwiseMover};
use bioshell_numerical::{get_random_point_nearby};
use bioshell_sim::{Energy, ResizableSystem, System};

use crate::{CartesianSystem};

pub struct RandomChain {
    pub bond_length: f64,
    pub energy_cutoff: f64,
    pub n_attempts: i16,
}

impl Default for RandomChain {
    fn default() -> Self {
        RandomChain {
            bond_length: 3.8,
            energy_cutoff: 0.00001,
            n_attempts: 100,
        }
    }
}

impl<E: Energy<CartesianSystem>> StepwiseMover<CartesianSystem, E> for RandomChain
{
    fn start(&mut self, system: &mut CartesianSystem, _energy: &E) -> f64
    {
        let c = system.get_box_len() / 2.0;
        system.set_size(2);
        system.set_xyz(0, c, c, c);
        let v = get_random_point_nearby(&system.get_coordinates()[0], self.bond_length);
        system.copy_from_vec(1, &v);

        return 1.0;
    }

    fn grow_by_one(&mut self, system: &mut CartesianSystem, energy: &E) -> f64 {
        let i = system.get_size();
        system.set_size(i + 1);
        let mut n_try = 0;
        while n_try < self.n_attempts {
            let v = get_random_point_nearby(&system.get_coordinates()[i - 1], self.bond_length);
            system.set_xyz(i, v.x, v.y, v.z);
            system.update_nbl(i);

            let en = energy.energy_by_pos(system, i);
            if en <= self.energy_cutoff {
                return 1.0;
            }
            n_try += 1;
        }
        return 0.0;
    }
}

impl<E: Energy<CartesianSystem>> StepwiseBuilder<CartesianSystem, E> for RandomChain {
    fn build(&mut self, system: &mut CartesianSystem, energy: &E) -> f64 {
        let mut step = RandomChain::default();
        step.start(system, energy);
        while system.get_size() < system.get_capacity() {
            step.grow_by_one(system, energy);
        }
        return 1.0;
    }
}