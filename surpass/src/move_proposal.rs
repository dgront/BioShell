use crate::SurpassAlphaSystem;

pub trait Mover<const N_MOVED: usize> {
    fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<N_MOVED>);
}

#[derive(Clone, Debug)]
pub struct MoveProposal<const N: usize> {
    pub first_moved_pos: usize,
    pub moved_cax: [i32; N],
    pub moved_cay: [i32; N],
    pub moved_caz: [i32; N],
}

impl<const N: usize> MoveProposal<N> {

    /// Creates new placeholder for proposed move coordinates, which are all set to 0.
    pub fn new() -> MoveProposal<N> {
        MoveProposal{
            first_moved_pos: 0,
            moved_cax: [0; N],
            moved_cay: [0; N],
            moved_caz: [0; N]
        }
    }

    pub fn apply(&self, model: &mut SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..N {
            model.cax[i_chain] = self.moved_cax[i_moved];
            model.cay[i_chain] = self.moved_cay[i_moved];
            model.caz[i_chain] = self.moved_caz[i_moved];
            i_chain += 1;
        }
    }

}
