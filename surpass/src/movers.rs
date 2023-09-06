use std::ops::Range;
use crate::SurpassAlphaSystem;

#[derive(Clone, Debug, Default)]
pub struct MoveProposal<const N: usize> {
    moved_range: Range<usize>,
    moved_cax: [i32; N],
    moved_cay: [i32; N],
    moved_caz: [i32; N],
}

impl<const N: usize> MoveProposal<N> {
    pub fn energy(&self) -> f64 { 0.0 }
}

pub struct HingeMove {
}

impl HingeMove {
    pub fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<8>) {}
}

pub struct IsothermalProtocol {
    temperature: f64,
    inner_steps: usize,
    outer_steps: usize,
}

impl IsothermalProtocol {
    pub fn run(&self, system: &mut SurpassAlphaSystem) {
        let mover = HingeMove{};
        let mut prp: MoveProposal<8> = MoveProposal::default();
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