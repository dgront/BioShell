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
        let mut rng = rand::thread_rng();
        let i_moved = rng.gen_range(0..system.size());

        let old_x: f64 = system.coordinates()[i_moved].x;
        let old_y: f64 = system.coordinates()[i_moved].y;
        let old_z: f64 = system.coordinates()[i_moved].z;
        let old_en: f64 = energy.energy_by_pos(system, i_moved);

        system.add(
            i_moved,
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
        );

        let new_en: f64 = energy.energy_by_pos(system, i_moved);
        if acc.check(old_en, new_en) {
            system.update_nbl(i_moved);
            self.succ_rate.n_succ += 1;
            return Option::from(i_moved..i_moved);
        } else {
            self.succ_rate.n_failed += 1;
            system.set(i_moved, old_x, old_y, old_z);
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

use bioshell_numerical::matrix::Matrix3x3;
use bioshell_numerical::vec3::Vec3;
use bioshell_numerical::Rototranslation;

/// A mover that applies a crankshaft move to a molecule.
pub struct CrankshaftMove {
    max_angle: f64,
    frag_size: usize,
    succ_rate: AcceptanceStatistics,
}

impl CrankshaftMove {
    /// Create a new mover that applies a crankshaft move.
    pub fn new(max_angle: f64, max_displacement: f64) -> CrankshaftMove {
        CrankshaftMove {
            max_angle,
            frag_size:5,
            succ_rate: Default::default(),
        }
    }

    fn apply(&self, system: &mut CartesianSystem, i1: usize, i2: usize, angle: f64) {
        let mut center = Vec3::zero();
        let mut end = Vec3::zero();
        center.add(&system.coordinates()[i1]);
        center.add(&system.coordinates()[i2]);
        center.mul(0.5);

        let v1 = system.coordinates()[i1].clone();
        let v2 = system.coordinates()[i2].clone();
        end.add(&v2);
        end.sub(&v1);

        let rt = Rototranslation::around_axis(&v1, &v2, angle);

        let mut v3 = rt.apply(&end);
        v3.sub(&end);
        v3.mul(self.max_displacement);

        let coords = system.coordinates();
        for i in i1+1..i2 {
                rt.apply_mut(&coords[i]);
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
        let size = system.size();

        let i1 = rng.gen_range(0..size-&self.frag_size-1);
        let mut i2 = i1 + &self.frag_size+1;

        let old_en = energy.energy(system);
        let angle = rng.gen_range(-&self.max_angle..&self.max_angle);
        self.apply(system, i1, i2, angle);
        let new_en = energy.energy(system);

        if acc.check(old_en, new_en) {
            self.succ_rate.n_succ += 1;
            return Option::from(i1..i2 + 1);
        } else {
            self.apply(system, i1, i2, -angle);
            return None;
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
