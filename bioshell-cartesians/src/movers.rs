use std::ops::Range;
use rand::Rng;

use crate::{CartesianSystem};
use bioshell_sim::System;
use bioshell_montecarlo::{AcceptanceStatistics, Mover};

/// performs a volume change
pub struct ChangeVolume {
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

impl ChangeVolume {
    pub fn new() -> ChangeVolume { ChangeVolume{ max_step: 0.1, succ_rate: Default::default() } }
}

impl Mover<CartesianSystem> for ChangeVolume {
    fn perturb(&mut self, system: &mut CartesianSystem) -> Range<usize> {
        let mut rng = rand::thread_rng();

        let v0 = system.volume();
        let lnV0 = v0.ln();
        let lnV = lnV0 + rng.gen_range(-self.max_step..self.max_step);
        let new_len = lnV.exp().powf(0.333333333);
        system.set_box_len(new_len);

        0..system.size()
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn add_success(&mut self) { self.succ_rate.n_succ += 1; }

    fn add_failure(&mut self) { self.succ_rate.n_failed += 1; }

    fn max_range(&self) -> f64 { self.max_step }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}

pub struct SingleAtomMove {
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

impl SingleAtomMove {
    pub fn new() -> SingleAtomMove { SingleAtomMove{ max_step: 2.0, succ_rate: Default::default() } }
}

impl Mover<CartesianSystem> for SingleAtomMove {
    fn perturb(&mut self, system: &mut CartesianSystem) -> Range<usize> {
        let mut rng = rand::thread_rng();
        let i_moved = rng.gen_range(0..system.size());
        system.add(i_moved,rng.gen_range(-self.max_step..self.max_step),
                   rng.gen_range(-self.max_step..self.max_step),rng.gen_range(-self.max_step..self.max_step));

        i_moved..i_moved
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn add_success(&mut self) { self.succ_rate.n_succ += 1; }

    fn add_failure(&mut self) { self.succ_rate.n_failed += 1; }

    fn max_range(&self) -> f64 { self.max_step }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}

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



