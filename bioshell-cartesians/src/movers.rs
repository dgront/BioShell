use std::ops::Range;
use rand::Rng;

use crate::{CartesianSystem};
use bioshell_sim::{Energy, System};
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};

/// A mover that moves a single, randomly selected atom by a small random vector.
///
/// Each move randomly selects an atom that is moved by `$\delta x,\delta y,\delta z$`.
/// The move is accepted or rejected according to any  [`AcceptanceCriterion`](AcceptanceCriterion)
/// struct defined in the **bioshell** library.
pub struct SingleAtomMove {
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

impl SingleAtomMove {
    /// Create a new mover that shifts a single atom.
    ///
    /// Each coordinate of  a randomly selected atom will be changed by at most `max_range`
    pub fn new(max_range: f64) -> SingleAtomMove { SingleAtomMove{ max_step: max_range, succ_rate: Default::default() } }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for SingleAtomMove {

    fn perturb(&mut self, system: &mut CartesianSystem, energy: &E, acc: &mut dyn AcceptanceCriterion) -> Option<Range<usize>> {

        let mut rng = rand::thread_rng();
        let i_moved = rng.gen_range(0..system.size());

        let old_x: f64 = system.coordinates()[i_moved].x;
        let old_y: f64 = system.coordinates()[i_moved].y;
        let old_z: f64 = system.coordinates()[i_moved].z;
        let old_en: f64 = energy.energy_by_pos(system, i_moved);

        system.add(i_moved,rng.gen_range(-self.max_step..self.max_step),
                   rng.gen_range(-self.max_step..self.max_step),rng.gen_range(-self.max_step..self.max_step));

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

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn max_range(&self) -> f64 { self.max_step }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}


/// A mover that changes a volume of a Cartesian system.
///
/// Each move changes volume of the simulation box by making a linear step in `$\ln(V)$` by `$\pm \delta$`.
/// A successful move changes the length of a cubic box and rescales accordingly atomic positions.
/// A move is accepted with probability `$\min(1, e^\omega)$`, where
///
/// ```math
/// \omega = \frac{\Delta E + p \Delta V}{kT} + (N+1)\ln(\frac{V_n}{V_n})
/// ```.
/// Such an implementation of the move leads to larger steps at larger volumes and smaller steps at smaller volumes
pub struct ChangeVolume {
    /// pressure of the system `$p$`
    pub pressure: f64,
    pub temperature: f64,
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

#[allow(non_upper_case_globals)]
const pV_to_Kelvins: f64 = 1.0E-30 / 1.380649E-23; // 10^-30 divided by the Boltzmann constant

impl ChangeVolume {

    pub fn new(pressure: f64, temperature: f64) -> ChangeVolume {
        ChangeVolume {pressure, temperature, max_step: 0.01, succ_rate: Default::default()}
    }
}


impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for ChangeVolume {

    #[allow(non_snake_case)]
    fn perturb(&mut self, system: &mut CartesianSystem, energy: &E, _acc: &mut dyn AcceptanceCriterion) -> Option<Range<usize>> {

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
        let weight = ((en_after - en_before) + self.pressure*(new_V - v0) * pV_to_Kelvins) / self.temperature -
            (system.size() + 1) as f64  * (new_V / v0).ln();

        if rng.gen_range(0.0..1.0) > (-weight/self.temperature).exp() { // --- move rejected
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


/*
pub struct PerturbChainFragment {
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

impl PerturbChainFragment {
    pub fn new() -> PerturbChainFragment { PerturbChainFragment{ max_step: 2.0, succ_rate: Default::default() } }
}

impl Mover<CartesianSystem> for PerturbChainFragment {

    fn perturb(&mut self, system: &mut CartesianSystem) -> Range<usize> {

        const N: usize = 3;
        const F: f64 = 2.0 / (1.0 + N as f64);

        let mut rng = rand::thread_rng();
        let mut moved_from = rng.gen_range(0..system.size()-N);
        let mut moved_to = moved_from + N - 1;
        while system.coordinates()[moved_from].chain_id != system.coordinates()[moved_to].chain_id {
            moved_from = rng.gen_range(0..system.size() - N);
            moved_to = moved_from + N;
        }

        let dx: f64 = rng.gen_range(-self.max_step..self.max_step) * F;
        let dy: f64 = rng.gen_range(-self.max_step..self.max_step) * F;
        let dz: f64 = rng.gen_range(-self.max_step..self.max_step) * F;

        for i in 0..N/2 {
            let fi: f64 = (i + 1) as f64;
            system.add(moved_from + i, dx * fi, dy * fi, dz * fi);
            system.add(moved_to - i, dx * fi, dy * fi, dz * fi);
        }

        if N % 2 == 1 {
            let mid = (moved_from + moved_to) / 2;
            system.add(mid,dx / F, dy / F, dz / F);
        }

        moved_from..moved_to
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn add_success(&mut self) { self.succ_rate.n_succ += 1; }

    fn add_failure(&mut self) { self.succ_rate.n_failed += 1; }

    fn max_range(&self) -> f64 { self.max_step }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}

*/

