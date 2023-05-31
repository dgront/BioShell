use rand::Rng;
use std::ops::Range;
use crate::CartesianSystem;
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Energy, System};
use bioshell_numerical::vec3::Vec3;
use bioshell_numerical::{random_point_nearby, Rototranslation};

/// A mover that changes a volume of a Cartesian system.
pub struct ChangeVolume {
    /// pressure of the system `$p$`
    pub pressure: f64,
    pub temperature: f64,
    max_step: f64,
    succ_rate: AcceptanceStatistics,
}

#[allow(non_upper_case_globals)]
const pV_to_Kelvins: f64 = 1.0E-30 / 1.380649E-23; // 10^-30 divided by the Boltzmann constant

impl ChangeVolume {
    pub fn new(pressure: f64, temperature: f64) -> ChangeVolume {
        ChangeVolume {
            pressure,
            temperature,
            max_step: 0.01,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for ChangeVolume {
    #[allow(non_snake_case)]
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        _acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
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
        let weight = ((en_after - en_before) + self.pressure * (new_V - v0) * pV_to_Kelvins)
            / self.temperature
            - (system.size() + 1) as f64 * (new_V / v0).ln();

        if rng.gen_range(0.0..1.0) > (-weight / self.temperature).exp() {
            // --- move rejected
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