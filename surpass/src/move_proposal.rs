use rand::Rng;
use crate::SurpassAlphaSystem;

pub trait Mover<const N_RES_MOVED: usize, const N_ATOMS_MOVED: usize> {
    fn propose<R: Rng>(&self, system: &SurpassAlphaSystem, rndgen: &mut R, proposal: &mut MoveProposal<N_RES_MOVED, N_ATOMS_MOVED>);
}

#[derive(Clone, Debug)]
pub struct MoveProposal<const N_RES: usize, const N_ATOMS: usize> {
    pub first_moved_pos: usize,
    pub bbx: [i32; N_ATOMS],
    pub bby: [i32; N_ATOMS],
    pub bbz: [i32; N_ATOMS],
}

impl<const N_RES: usize, const N_ATOMS: usize> MoveProposal<N_RES, N_ATOMS> {

    /// Creates new placeholder for proposed move coordinates, which are all set to 0.
    pub fn new() -> MoveProposal<N_RES, N_ATOMS> {
        MoveProposal{
            first_moved_pos: 0,
            bbx: [0; N_ATOMS],
            bby: [0; N_ATOMS],
            bbz: [0; N_ATOMS]
        }
    }

    pub fn apply(&self, model: &mut SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..N_ATOMS {
            model.bbx[i_chain] = self.bbx[i_moved];
            model.bby[i_chain] = self.bby[i_moved];
            model.bbz[i_chain] = self.bbz[i_moved];
            i_chain += 1;
        }
    }

    pub fn backup(&mut self, model: &SurpassAlphaSystem) {
        let mut i_chain = self.first_moved_pos;
        for i_moved in 0..N_ATOMS {
            self.bbx[i_moved] = model.bbx[i_chain];
            self.bby[i_moved] = model.bby[i_chain];
            self.bbz[i_moved] = model.bbz[i_chain];
            i_chain += 1;
        }
    }
}
