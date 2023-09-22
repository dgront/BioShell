use crate::{HingeMove, MoveProposal, Mover, SurpassAlphaSystem};

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