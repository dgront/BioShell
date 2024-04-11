use rand::Rng;
use crate::SurpassAlphaSystem;

const N_MOVED_MAX: usize = 8;

pub trait Mover {
    fn propose<R: Rng>(&self, system: &SurpassAlphaSystem, rndgen: &mut R, proposal: &mut MoveProposal);
}

#[derive(Clone, Debug)]
pub struct MoveProposal {
    pub first_moved_pos: usize,
    pub n_moved: usize,
    pub cax: [i32; N_MOVED_MAX],
    pub cay: [i32; N_MOVED_MAX],
    pub caz: [i32; N_MOVED_MAX],
}

impl MoveProposal {

    /// Creates new placeholder for proposed move coordinates, which are all set to 0.
    pub fn new(n_moved: usize) -> MoveProposal {
        MoveProposal{
            first_moved_pos: 0, n_moved,
            cax: [0; N_MOVED_MAX], cay: [0; N_MOVED_MAX], caz: [0; N_MOVED_MAX]
        }
    }

    /// Alters the given model by copying into it the coordinates from this move proposal
    pub fn apply(&self, model: &mut SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..self.n_moved {
            model.cax[i_chain] = self.cax[i_moved];
            model.cay[i_chain] = self.cay[i_moved];
            model.caz[i_chain] = self.caz[i_moved];
            i_chain += 1;
        }
    }

    /// Alters this move proposal by copying coordinates from the given model
    pub fn backup(&mut self, model: &SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..self.n_moved {
            self.cax[i_moved] = model.cax[i_chain];
            self.cay[i_moved] = model.cay[i_chain];
            self.caz[i_moved] = model.caz[i_chain];
            i_chain += 1;
        }
    }
}
