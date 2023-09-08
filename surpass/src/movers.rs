use std::ops::Range;
use rand::{Rng, thread_rng};
use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::SurpassAlphaSystem;

#[derive(Clone, Debug)]
pub struct MoveProposal<const N: usize> {
    moved_range: Range<usize>,
    moved_cax: [i32; N],
    moved_cay: [i32; N],
    moved_caz: [i32; N],
}

impl<const N: usize> MoveProposal<N> {

    /// Creates new placeholder for proposed move coordinates, which are all set to 0.
    pub fn new() -> MoveProposal<N> {
        MoveProposal{
            moved_range: Default::default(),
            moved_cax: [0; N],
            moved_cay: [0; N],
            moved_caz: [0; N]
        }
    }

    pub fn energy(&self) -> f64 { 0.0 }
}

const HINGE_MOVE_SIZE: usize = 8;

pub struct HingeMove {
    max_angle: f64,
    max_range_allowed: f64
}

impl HingeMove {
    pub fn new(max_angle: f64, max_range_allowed: f64) -> HingeMove {
        HingeMove{ max_angle,max_range_allowed }
    }

    pub fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<HINGE_MOVE_SIZE>) {
        let mut rng = thread_rng();
        // --- pick end points randomly
        let mut moved_from = 0;
        let mut moved_to = HINGE_MOVE_SIZE + 1;
        loop {
            let moved_from = rng.gen_range(0..system.count_atoms() - HINGE_MOVE_SIZE - 1);
            let moved_to = moved_from + HINGE_MOVE_SIZE - 1;
            if system.chain(moved_from) == system.chain(moved_to) { break };
        }
        // --- prepare rototranslation
        let start_vector: Vec3 = Vec3::new(system.cax(moved_from), system.cay(moved_from), system.caz(moved_from));
        let end_vector: Vec3 = Vec3::new(system.cax(moved_to), system.cay(moved_to), system.caz(moved_to));
        let angle = rng.gen_range(-self.max_angle..self.max_angle);
        let roto = Rototranslation::around_axis(&start_vector, &end_vector, angle);
        // --- rotate atoms around the axis
    }
}

pub struct IsothermalProtocol {
    temperature: f64,
    inner_steps: usize,
    outer_steps: usize,
}

impl IsothermalProtocol {
    pub fn run(&self, system: &mut SurpassAlphaSystem) {
        let mover = HingeMove::new(30.0_f64.to_radians(), 60.0_f64.to_radians());
        let mut prp: MoveProposal<8> = MoveProposal::new();
        for o in 0..self.outer_steps {
            for i in 0..self.inner_steps {
                // --- propose a move
                mover.propose(&system, &mut prp);
                // --- evaluate deltaE
                // --- check Metropolis criterion
                // --- apply the move if successful (40% chance)
            }
        }
    }
}