use std::ops::Range;

use crate::ff::Energy;
use crate::{System, ZeroEnergy};

/// A pairwise nonbonded interaction can evaluate its energy just from a squared distance value
/// between two atoms
pub trait PairwiseNonbonded {

    fn energy_for_distance_squared(&self, d2: f32) -> f64;
}

macro_rules! pairwise_contact_kernel {
    ($x:expr, $y:expr, $z:expr, $chain:expr, $i:expr, $self:expr, $en:expr) => {

        let mut d = $chain.delta_x($i, $x);
        let mut d2 = d*d;
        d = $chain.delta_y($i, $y);
        d2 += d*d;
        d = $chain.delta_z($i, $z);
        d2 += d*d;
        $en += $self.energy.energy_for_distance_squared(d2);
    }
}

macro_rules! pairwise_contact_neighbors_loop {
    ($self:expr, $chain:expr, $pos:expr, $neighbors:expr, $en:expr) => {
            let x: f32 = $chain[$pos].x;
            let y: f32 = $chain[$pos].y;
            let z: f32 = $chain[$pos].z;
            for i in $neighbors {
                pairwise_contact_kernel!(x, y, z, $chain, *i, $self, $en);
            }
    }
}


/// Evaluates pairwise nonbonded interactions between all atoms.
///
/// The evaluator uses the provided [PairwiseNonbonded] instance to evaluate energy between all atoms
/// in a given system.
#[allow(unused)]
pub struct PairwiseNonbondedEvaluator {

    pub cutoff: f32,
    cutoff_square: f32,
    pub energy: Box<dyn PairwiseNonbonded>,
}

impl PairwiseNonbondedEvaluator {

    pub fn new(cutoff: f32, pairwise_energy: Box<dyn PairwiseNonbonded>) -> PairwiseNonbondedEvaluator {
        PairwiseNonbondedEvaluator { cutoff: cutoff, cutoff_square:cutoff*cutoff, energy: pairwise_energy }
    }
}

impl Energy for PairwiseNonbondedEvaluator {

    fn energy(&self, system: &System) -> f64 {
        let mut en:f64 = 0.0;
        for i in 0..system.size() {
            en += self.energy_by_pos(system, i);
        }

        return en / 2.0;
    }

    fn energy_by_pos(&self, system: &System, pos:usize) -> f64 {
        let mut en: f64 = 0.0;

        let chain = system.coordinates();
        let pos_neighbors = system.neighbor_list().neighbors(pos);
        pairwise_contact_neighbors_loop!(self, chain, pos, pos_neighbors, en);

        return en;
    }

    fn delta_energy_by_range(&self, old: &System, moved: &Range<usize>, new: &System) -> f64 {

        let mut en_old: f64 = 0.0;
        let mut en_new: f64 = 0.0;
        let old_coords = old.coordinates();
        let new_coords = new.coordinates();
        for jpos in moved.start..moved.end + 1 {
            // --- compute energy before and after a move for each moved jpos position
            pairwise_contact_neighbors_loop!(self, old_coords, jpos, old.neighbor_list().neighbors(jpos), en_old);
            pairwise_contact_neighbors_loop!(self, new_coords, jpos, new.neighbor_list().neighbors(jpos), en_new);
        }

        return en_new - en_old;
    }

    fn name(&self) -> String { String::from("PairwiseNonbonded") }
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