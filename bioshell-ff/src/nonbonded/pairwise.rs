use std::ops::Range;

use bioshell_sim::{System, Energy};
use bioshell_cartesians::{CartesianSystem};

/// A pairwise nonbonded interaction can evaluate its energy just from a squared distance value
/// between two atoms
pub trait PairwiseNonbonded {

    fn energy_for_distance_squared(&self, d2: f64) -> f64;
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
            let x: f64 = $chain[$pos].x;
            let y: f64 = $chain[$pos].y;
            let z: f64 = $chain[$pos].z;
            for i in $neighbors {
                pairwise_contact_kernel!(x, y, z, $chain, *i, $self, $en);
            }
    }
}


/// Evaluates pairwise nonbonded interactions between all atoms.
///
/// The evaluator uses the provided [PairwiseNonbonded] instance to evaluate energy between all atoms
/// in a given system.
pub struct PairwiseNonbondedEvaluator<E:PairwiseNonbonded> {

    pub cutoff: f64,
    cutoff_square: f64,
    pub energy: E,
}

impl<E:PairwiseNonbonded> PairwiseNonbondedEvaluator<E> {

    pub fn new(cutoff: f64, pairwise_energy: E) -> PairwiseNonbondedEvaluator<E> {
        PairwiseNonbondedEvaluator { cutoff: cutoff, cutoff_square:cutoff*cutoff, energy: pairwise_energy }
    }
}

impl<E:PairwiseNonbonded> Energy<CartesianSystem> for PairwiseNonbondedEvaluator<E> {

    fn energy(&self, system: &CartesianSystem) -> f64 {
        let mut en:f64 = 0.0;
        for i in 0..system.size() {
            en += self.energy_by_pos(system, i);
        }

        return en / 2.0;
    }

    fn energy_by_pos(&self, system: &CartesianSystem, pos:usize) -> f64 {
        let mut en: f64 = 0.0;

        let chain = system.coordinates();
        let pos_neighbors = system.neighbor_list().neighbors(pos);
        pairwise_contact_neighbors_loop!(self, chain, pos, pos_neighbors, en);

        return en;
    }

    fn energy_by_range(&self, system: &CartesianSystem, range: &Range<usize>) -> f64 {
        todo!()
    }


    fn name(&self) -> String { String::from("PairwiseNonbonded") }
}
//----------------


pub struct SimpleContact {
    /// repulsion distance
    pub repulsion_distance: f64,
    /// energy well starts
    pub contact_from_distance: f64,
    /// energy well end
    pub contact_to_distance: f64,
    /// energy value
    pub contact_energy: f64,
    /// repulsion value
    pub repulsion_energy: f64,

    r_c_sq: f64,
    r_from_sq: f64,
    r_to_sq: f64,
}


impl SimpleContact {
    pub fn new(r_repulsion: f64, r_from: f64, r_to: f64, en_repulsion: f64, en_contact: f64) -> SimpleContact {
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
    fn energy_for_distance_squared(&self, d2: f64) -> f64 {
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

    fn energy_for_distance_squared(&self, r2: f64) -> f64 {
        return if (r2 as f64) < self.cutoff_sq {
            let r2_s: f64 = self.sigma_sq / r2 as f64;
            let r6: f64 = r2_s * r2_s * r2_s;
            let r12: f64 = r6 * r6;
            self.epsilon_four * (r12 - r6)
        } else { 0.0 }
    }
}