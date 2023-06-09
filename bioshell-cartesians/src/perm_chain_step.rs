//! Provides means to build a Cartesian system, either in a deterministic or a stochastic way
use rand::Rng;

use bioshell_montecarlo::{StepwiseMover};
use bioshell_numerical::{get_random_point_nearby, Vec3};
use bioshell_sim::{Energy, ResizableSystem, System};

use crate::{Coordinates};



pub struct PERMChainStep {
    pub temperature: f64,
    pub bond_length: f64,
    pub n_trials: i16,
}

impl PERMChainStep {
    pub fn new(temperature: f64, n_trials: i16) -> PERMChainStep {
        PERMChainStep {
            temperature,
            bond_length: 3.8,
            n_trials,
        }
    }
}

impl<E: Energy<Coordinates>> StepwiseMover<Coordinates, E> for PERMChainStep {
    /// Starts a new chain by placing its first bead in the center of a simulation box
    ///
    /// Always returns 1.0 for the statistical weight of the newly started chain as a single bead
    /// has nothing to interact with. The energy parameter is not used therefore.
    fn start(&mut self, system: &mut Coordinates, _energy: &E) -> f64 {
        let c = system.get_box_len() / 2.0;
        system.set_xyz(0, c, c, c);
        system.set_size(1);

        return 1.0;
    }

    fn grow_by_one(&mut self, system: &mut Coordinates, energy: &E) -> f64 {
        let i = system.get_size();
        system.set_size(i + 1);

        let mut weights: Vec<f64> = Vec::with_capacity(self.n_trials as usize);
        let mut vn: Vec<Vec3> = Vec::with_capacity(self.n_trials as usize);
        let center: Vec3 = (&system[i - 1]).clone();
        // ---------- propose n_trials random proposals and score them
        for _k in 0..self.n_trials {
            let v_k = get_random_point_nearby(&center, self.bond_length);
            system.from_vec(i, &v_k);
            let en = energy.energy_by_pos(system, i);
            weights.push((-en / self.temperature).exp());
            vn.push(v_k);
        }
        let total = weights.iter().sum();
        if total < 1e-100 {
            return 0.0;
        } // --- no suitable move generated

        // ---------- select one of the possible extension by importance sampling
        let mut rng = rand::thread_rng();
        let r = rng.gen_range(0.0..total);
        let mut which_v: usize = 0;
        let mut s = weights[which_v];
        while s <= r {
            which_v += 1;
            s += weights[which_v]
        }

        // ---------- set the coordinates
        system.set_xyz(i, vn[which_v].x, vn[which_v].y, vn[which_v].z);

        // ---------- return the statistical weight
        return total;
    }
}