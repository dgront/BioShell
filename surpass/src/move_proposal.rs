use crate::SurpassAlphaSystem;

pub trait Mover<const N_MOVED: usize> {
    fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<N_MOVED>);
}

#[derive(Clone, Debug)]
pub struct MoveProposal<const N: usize> {
    pub first_moved_pos: usize,
    pub cax: [i32; N],
    pub cay: [i32; N],
    pub caz: [i32; N],
}

impl<const N: usize> MoveProposal<N> {

    /// Creates new placeholder for proposed move coordinates, which are all set to 0.
    pub fn new() -> MoveProposal<N> {
        MoveProposal{
            first_moved_pos: 0,
            cax: [0; N],
            cay: [0; N],
            caz: [0; N]
        }
    }

    pub fn apply(&self, model: &mut SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..N {
            model.cax[i_chain] = self.cax[i_moved];
            model.cay[i_chain] = self.cay[i_moved];
            model.caz[i_chain] = self.caz[i_moved];
            i_chain += 1;
        }
    }

    pub fn backup(&mut self, model: &SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..N {
            self.cax[i_moved] = model.cax[i_chain];
            self.cay[i_moved] = model.cay[i_chain];
            self.caz[i_moved] = model.caz[i_chain];
            i_chain += 1;
        }
    }
}
