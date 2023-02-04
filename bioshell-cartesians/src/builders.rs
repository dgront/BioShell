use bioshell_montecarlo::{StepwiseBuilder, StepwiseMover};

use crate::{CartesianSystem, Coordinates};
use bioshell_sim::{System, Energy};
use bioshell_numerical::{random_point_nearby, random_unit_versor, Vec3};

pub struct RandomChain {
    pub bond_length: f64,
    pub energy_cutoff: f64,
    pub n_attempts: i16
}

impl Default for RandomChain {
    fn default() -> Self { RandomChain {bond_length: 3.8, energy_cutoff:0.00001, n_attempts: 100} }
}

impl<E: Energy<CartesianSystem>> StepwiseMover<CartesianSystem, E> for RandomChain {

    fn start(&mut self, system: &mut CartesianSystem, energy: &E) -> f64 {
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

