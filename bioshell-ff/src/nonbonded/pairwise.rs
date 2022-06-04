use std::ops::Range;

use crate::ff::Energy;
use crate::{Coordinates, ZeroEnergy};

/// A pairwise nonbonded interaction cat evaluate its energy just from a squared distance value
pub trait PairwiseNonbonded {

    fn energy_for_distance_squared(&self, d2: f32) -> f64;
}

macro_rules! pairwise_contact_kernel {
    ($x:expr, $y:expr, $z:expr, $chain:expr, $i:expr, $self:expr, $en:expr) => {

        let mut d = $chain.delta_x($i, $x);
        let mut d2 = d*d;
        if d2 > $self.cutoff_square { continue; }
        d = $chain.delta_y($i, $y);
        d2 += d*d;
        if d2 > $self.cutoff_square { continue; }
        d = $chain.delta_z($i, $z);
        d2 += d*d;
        $en += $self.energy.energy_for_distance_squared(d2);
    }
}

macro_rules! pairwise_contact_kernel_loop {
    ($start:expr, $end:expr, $x:expr, $y:expr, $z:expr, $chain:expr, $self:expr, $en:expr) => {
            for i in $start..$end {
                pairwise_contact_kernel!($x, $y, $z, $chain, i, $self, $en);
            }
    }
}


/// Evaluates pairwise nonbonded interactions
///
///
pub struct PairwiseNonbondedEvaluator {

    /// sequence separation, e.g. 0 for argon fluid, 1 for simple polymer
    pub sequence_separation: usize,
    pub cutoff: f32,
    cutoff_square: f32,
    pub energy: Box<dyn PairwiseNonbonded>,
}

impl PairwiseNonbondedEvaluator {

    pub fn new(separation: usize, cutoff: f32, pairwise_energy: Box<dyn PairwiseNonbonded>) -> PairwiseNonbondedEvaluator {
        PairwiseNonbondedEvaluator { sequence_separation: separation, cutoff: cutoff,
            cutoff_square:cutoff*cutoff, energy: pairwise_energy }
    }

    pub fn pair_energy(&self, system: &Coordinates, ipos: usize, jpos: usize) -> f64 {
        if ipos > jpos && ipos - jpos <= self.sequence_separation { return 0.0; }
        if ipos < jpos && jpos - ipos <= self.sequence_separation { return 0.0; }
        let d = system.closest_distance_square(ipos, jpos).sqrt();

        self.energy.energy_for_distance_squared(d)
    }

    fn each_vs_each_energy(&self, system: &Coordinates, moved: &Range<usize>) ->f64 {

        let mut en: f64 = 0.0;
        for ipos in moved.start..moved.end {
            let xi: f32 = system[ipos].x;
            let yi: f32 = system[ipos].y;
            let zi: f32 = system[ipos].z;

            for jpos in (ipos+1+self.sequence_separation)..moved.end + 1 {
                pairwise_contact_kernel!(xi, yi, zi, system, jpos, self, en);
            }
        }
        return en;
    }
}

impl Energy for PairwiseNonbondedEvaluator {

    fn energy(&self, system: &Coordinates) -> f64 {
        let mut en:f64 = 0.0;
        for i in 0..system.size() {
            en += self.energy_by_pos(system, i);
        }

        return en / 2.0;
    }

    fn energy_by_pos(&self, chain: &Coordinates, pos:usize) -> f64 {
        let mut en: f64 = 0.0;

        let x: f32 = chain[pos].x;
        let y: f32 = chain[pos].y;
        let z: f32 = chain[pos].z;

        if pos > self.sequence_separation {
            pairwise_contact_kernel_loop!(0, pos - self.sequence_separation, x, y, z, chain, self, en);
        }
        // for the chain of 10 atoms and sep=1, we score when pos is at most 7 (7 vs 9)
        if pos < chain.size() - 1 - self.sequence_separation {
            pairwise_contact_kernel_loop!(pos + 1 + self.sequence_separation, chain.size(), x, y, z, chain, self, en);
        }

        return en;
    }

    fn delta_energy_by_range(&self, old: &Coordinates, moved: &Range<usize>, new: &Coordinates) -> f64 {

        let mut en_old: f64 = 0.0;
        let mut en_new: f64 = 0.0;

        for jpos in moved.start..moved.end + 1 {
            // --- x, y, z of a moved atom are different between new and old systems
            let xo: f32 = old[jpos].x;
            let yo: f32 = old[jpos].y;
            let zo: f32 = old[jpos].z;
            let xn: f32 = new[jpos].x;
            let yn: f32 = new[jpos].y;
            let zn: f32 = new[jpos].z;

            // --- Energy upstream the moved range; e.g. for sep=1 we calculate energy here if moved 2 or greater
            if moved.start > self.sequence_separation {
                pairwise_contact_kernel_loop!(0, moved.start - self.sequence_separation, xo, yo, zo, old, self, en_old);
                pairwise_contact_kernel_loop!(0, moved.start - self.sequence_separation, xn, yn, zn, new, self, en_new);
            }
            // --- Energy downstream the moved range, e.g. for N=10 and sep=1 we calculate here if moved at most 7 (7 vs 9)
            if moved.end < old.size() - self.sequence_separation -1 {
                let start = moved.end + 1 + self.sequence_separation;
                pairwise_contact_kernel_loop!(start, old.size(), xo, yo, zo, old, self, en_old);
                pairwise_contact_kernel_loop!(start, old.size(), xn, yn, zn, new, self, en_new);
            }
        }
        // --- Energy within the moved range
        en_old += self.each_vs_each_energy(old,moved);
        en_new += self.each_vs_each_energy(new,moved);

        return en_new - en_old;
    }
}

impl PairwiseNonbonded for ZeroEnergy {
    fn energy_for_distance_squared(&self, _d2: f32) -> f64 { 0.0 }
}

pub struct SimpleContact {
    /// repulsion distance
    pub repulsion_distance: f32,
    /// energy well starts
    pub contact_from_distance: f32,
    /// energy well end
    pub contact_to_distance: f32,
    /// energy value
    pub contact_energy: f32,
    /// repulsion value
    pub repulsion_energy: f32,

    r_c_sq: f32,
    r_from_sq: f32,
    r_to_sq: f32,
}


impl SimpleContact {
    pub fn new(r_repulsion: f32, r_from: f32, r_to: f32, en_repulsion: f32, en_contact: f32) -> SimpleContact {
        SimpleContact {
            repulsion_distance: r_repulsion,
            contact_from_distance: r_from,
            contact_to_distance: r_to,
            repulsion_energy: en_repulsion,
            contact_energy: en_contact,
            r_c_sq: r_repulsion * r_repulsion,
            r_from_sq: r_from * r_from,
            r_to_sq: r_to * r_to
        }
    }
}

impl PairwiseNonbonded for SimpleContact {
    fn energy_for_distance_squared(&self, d2: f32) -> f64 {
        if d2 > self.r_to_sq { return 0.0; }
        if d2 < self.r_c_sq { return self.repulsion_energy as f64; }
        if d2 > self.r_from_sq { return self.contact_energy as f64; } else { return 0.0; }
    }
}

pub struct LennardJonesHomogenic {
    /// Energy constant
    pub epsilon: f64,
    /// VdW radius
    pub sigma: f64,
    /// Interaction cutoff
    pub cutoff: f64,
    sigma_sq: f64,
    cutoff_sq: f64,
    epsilon_four: f64,
}

impl LennardJonesHomogenic {
    pub fn new(epsilon: f64, sigma: f64, cutoff: f64 ) -> LennardJonesHomogenic {
            LennardJonesHomogenic {epsilon, sigma, cutoff,
                sigma_sq: sigma*sigma, cutoff_sq: cutoff*cutoff, epsilon_four: 4.0 * epsilon }
    }
}

impl PairwiseNonbonded for LennardJonesHomogenic {

    fn energy_for_distance_squared(&self, r2: f32) -> f64 {
        return if (r2 as f64) < self.cutoff_sq {
            let r2_s: f64 = self.sigma_sq / r2 as f64;
            let r6: f64 = r2_s * r2_s * r2_s;
            let r12: f64 = r6 * r6;
            self.epsilon_four * (r12 - r6)
        } else { 0.0 }
    }
}