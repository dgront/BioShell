use std::ops::Range;
use rand::Rng;

use crate::{CartesianSystem};
use bioshell_sim::{Energy, System};
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, MetropolisCriterion, Mover, Sampler};

/// Adds a volume changing move to a given protocol
pub struct VolumeChangingProtocol<E: Energy<CartesianSystem>> {
    pub pressure: f64,
    max_step: f64,
    succ_rate: AcceptanceStatistics,
    sampler: Box<dyn Sampler<CartesianSystem, E, MetropolisCriterion>>,
}


impl<E: Energy<CartesianSystem>> VolumeChangingProtocol<E> {

    pub fn new(pressure: f64, mut sampler: Box<dyn Sampler<CartesianSystem, E, MetropolisCriterion>>) -> VolumeChangingProtocol<E> {

        VolumeChangingProtocol {pressure, max_step: 0.01, succ_rate: Default::default(), sampler}
    }
}

impl< E: Energy<CartesianSystem>> Sampler<CartesianSystem, E, MetropolisCriterion>  for VolumeChangingProtocol<E> {

    fn make_sweeps(&mut self, n: usize, coords: &mut CartesianSystem, energy: &E) {

        // ---------- perform "regular" sampling
        self.sampler.make_sweeps(n, coords, energy);

    }

    fn acceptance_criterion(&mut self) -> &mut MetropolisCriterion { self.sampler.acceptance_criterion() }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E, MetropolisCriterion> for VolumeChangingProtocol<E> {

    fn perturb(&mut self, system: &mut CartesianSystem, energy: &E, acc: &mut MetropolisCriterion) -> Option<Range<usize>> {

        // ---------- attempt the volume change
        let en_before = energy.energy(system);
        let mut rng = rand::thread_rng();
        let v0 = system.volume();
        let lnV0 = v0.ln();
        let lnV = lnV0 + rng.gen_range(-self.max_step..self.max_step);
        let new_V = lnV.exp();
        let old_len = system.box_len();
        let new_len = new_V.powf(0.333333333);
        let f = new_len / system.box_len();
        for i in 0..system.size() {
            let x = system.coordinates().x(i) * f;
            let y = system.coordinates().y(i) * f;
            let z = system.coordinates().z(i) * f;
            system.set(i, x, y, z);
        }
        system.set_box_len(new_len);
        let en_after = energy.energy(system);
        let weight = (en_after - en_before) + self.pressure*(new_V - v0) -
            (system.size() + 1) as f64 * acc.temperature * (new_V / v0).ln();

        if rng.gen_range(0.0..1.0) > (-weight/acc.temperature).exp() { // --- move rejected
            for i in 0..system.size() {
                let x = system.coordinates().x(i) / f;
                let y = system.coordinates().y(i) / f;
                let z = system.coordinates().z(i) / f;
                system.set(i, x, y, z);
            }
            self.succ_rate.n_failed += 1;
            system.set_box_len(old_len);
            return None;
        } else {
            self.succ_rate.n_succ += 1;
            return Some(0..system.size());
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn max_range(&self) -> f64 { self.max_step }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}