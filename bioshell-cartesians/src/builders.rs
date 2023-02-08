use rand::{Rng};

use bioshell_montecarlo::{StepwiseBuilder, StepwiseMover};
use bioshell_sim::{System, Energy};
use bioshell_numerical::{random_point_nearby, Vec3};

use crate::{CartesianSystem, Coordinates};

pub struct RandomChain {
    pub bond_length: f64,
    pub energy_cutoff: f64,
    pub n_attempts: i16
}

impl Default for RandomChain {
    fn default() -> Self { RandomChain {bond_length: 3.8, energy_cutoff:0.00001, n_attempts: 100} }
}

impl<E: Energy<CartesianSystem>> StepwiseMover<CartesianSystem, E> for RandomChain {

    fn start(&mut self, system: &mut CartesianSystem, _energy: &E) -> f64 {
        let c = system.box_len() / 2.0;
        system.set(0, c, c, c);
        let v = random_point_nearby(&system.coordinates()[0], self.bond_length);
        system.copy_from_vec(1, &v);
        system.set_size(2);

        return 1.0
    }

    fn grow_by_one(&mut self, system: &mut CartesianSystem, energy: &E) -> f64 {

        let i = system.size();
        system.set_size(i + 1);
        let mut n_try = 0;
        while n_try < self.n_attempts {
            let v = random_point_nearby(&system.coordinates()[i-1], self.bond_length);
            system.copy_from_vec(i, &v);
            system.update_nbl(i);

            let en = energy.energy_by_pos(system, i);
            if en <= self.energy_cutoff {
                return 1.0;
            }
            n_try += 1;
        }
        return 0.0
    }
}

impl<E: Energy<CartesianSystem>> StepwiseBuilder<CartesianSystem, E> for RandomChain {

    fn build(&mut self, system: &mut CartesianSystem, energy: &E) -> f64 {
        let mut step = RandomChain::default();
        let mut w_total = 1.0;
        step.start(system, energy);
        while system.size() < system.capacity() {
            let w = step.grow_by_one(system, energy);
            w_total *= w;
        }
        return w_total;
    }
}


pub struct PERMChainStep {
    pub temperature: f64,
    pub bond_length: f64,
    pub n_trials: i16,
}

impl PERMChainStep {
    pub fn new(temperature: f64, n_trials: i16) -> PERMChainStep {
        PERMChainStep { temperature, bond_length: 3.8, n_trials }
    }
}

impl<E: Energy<Coordinates>> StepwiseMover<Coordinates, E> for PERMChainStep {

    /// Starts a new chain by placing its first bead in the center of a simulation box
    ///
    /// Always returns 1.0 for the statistical weight of the newly started chain as a single bead
    /// has nothing to interact with. The energy parameter is not used therefore.
    fn start(&mut self, system: &mut Coordinates, _energy: &E) -> f64 {
        let c = system.box_len() / 2.0;
        system.set(0, c, c, c);
        system.set_size(1);

        return 1.0
    }

    fn grow_by_one(&mut self, system: &mut Coordinates, energy: &E) -> f64 {

        let i = system.size();
        system.set_size(i + 1);

        let mut weights: Vec<f64> = Vec::with_capacity(self.n_trials as usize);
        let mut vn: Vec<Vec3> = Vec::with_capacity(self.n_trials as usize);
        let center: Vec3 = (&system[i-1]).clone();
        // ---------- propose n_trials random proposals and score them
        for k in 0..self.n_trials {
            let v_k = random_point_nearby(&center, self.bond_length);
            system.copy_from_vec(i, &v_k);
            let en = energy.energy_by_pos(system, i);
            weights.push((-en/self.temperature).exp());
            vn.push(v_k);
        }
        let total = weights.iter().sum();
        if total < 1e-100 { return 0.0; }       // --- no suitable move generated

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
        system.copy_from_vec(i, &vn[which_v]);

        // ---------- return the statistical weight
        return total;
    }
}

pub fn cubic_grid_atoms(system: &mut Coordinates) {

    let points_one_side: usize = (f64::powf(system.size() as f64, 1.0 / 3.0)).ceil() as usize;
    let dw = system.box_len() / points_one_side as f64;
    let cell_margin = dw / 2.0;

    for i in 0..system.size() {
        let k = i % points_one_side;
        let l = (i / points_one_side) % points_one_side;
        let m = (i / (points_one_side * points_one_side)) % points_one_side;
        system.set(i,dw * k as f64 + cell_margin,dw * l as f64 + cell_margin,dw * m as f64 + cell_margin)
    }
}


pub fn square_grid_atoms(system: &mut Coordinates) {

    let points_one_side: usize = (f64::powf(system.size() as f64, 0.5)).ceil() as usize;
    let dw = system.box_len() / points_one_side as f64;
    let cell_margin = dw / 2.0;

    for i in 0..system.size() {
        let k = i % points_one_side;
        let l = i / points_one_side;
        system.set(i,dw * k as f64 + cell_margin,dw * l as f64 + cell_margin,0.0);
    }
}

