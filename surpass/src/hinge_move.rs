#[warn(unused_imports)]
use log::{trace};
use rand::{Rng};
use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::{MoveProposal, Mover, SurpassAlphaSystem};

pub struct HingeMove<const N_RES_MOVED: usize, const N_ATOMS_MOVED: usize> {
    max_angle: f64,
    max_range_allowed: f64
}

impl<const N_RES_MOVED: usize, const N_ATOMS_MOVED: usize> Mover<N_RES_MOVED, N_ATOMS_MOVED> for HingeMove<N_RES_MOVED, N_ATOMS_MOVED> {
    fn propose<R: Rng>(&self, system: &SurpassAlphaSystem, rnd_gen: &mut R, proposal: &mut MoveProposal<N_RES_MOVED, N_ATOMS_MOVED>) {

        // --- pick end points randomly
        let mut moved_from;
        let mut moved_to;
        loop {
            moved_from = rnd_gen.gen_range(1..system.count_residues() - N_RES_MOVED);
            moved_to = moved_from + N_RES_MOVED - 1;
            if system.chain(moved_from-1) == system.chain(moved_to+1) { break };
        }
        let angle = rnd_gen.gen_range(-self.max_angle..self.max_angle);
        trace!("hinge move of {}:{} by {}", moved_from, moved_from+N_RES_MOVED-1, angle);
        self.compute_move(system, moved_from, angle, proposal);
    }
}

impl<const N_RES_MOVED: usize, const N_ATOMS_MOVED: usize> HingeMove<N_RES_MOVED, N_ATOMS_MOVED> {

    pub fn new(max_angle: f64, max_range_allowed: f64) -> HingeMove<N_RES_MOVED, N_ATOMS_MOVED> {
        HingeMove{ max_angle, max_range_allowed }
    }


    /// Computes moved coordinates by applying the appropriate rotation.
    ///
    /// # Arguments
    ///  * `system` - system to be modified
    ///  * `moved_from` - index of the first moved atom; must be greater than 1 and smaller than `N-8` where
    ///     `N` is the total number of atoms in this `system`. **Note** that this method doesn't check
    ///     whether the fragment subjected to move is located within a single chain.
    ///  * `angle` - rotation angle
    ///  * `proposal` - object used to store the resulting (rotated) coordinates
    pub fn compute_move(&self, system: &SurpassAlphaSystem, moved_from: usize, angle: f64, proposal: &mut MoveProposal<N_RES_MOVED, N_ATOMS_MOVED>) {

        // ------ Working example ----
        // Moving 8 atoms indexed from 1 to 8 right-exclusive (i.e. moved_from = 1) assumes:
        //  - whole window is 10 atoms long
        //  - pivot positions are 0 and 9
        //  - rotated residue range is 1..9 (9th is not rotated)
        //  - rotated atoms: 1*4..9*4
        // -------------------------
        // --- prepare rototranslation
        let start_vector: Vec3 = system.atom_to_vec3(moved_from - 1);
        let end_vector: Vec3 = system.atom_to_nearest_vec3(moved_from + N_RES_MOVED, moved_from - 1);
        let roto = Rototranslation::around_axis(&start_vector, &end_vector, angle);
        // --- rotate atoms around the axis and copy outcome to a proposal
        proposal.first_moved_pos = moved_from;
        let mut v: Vec3;
        let mut i_chain = moved_from*4;
        let mut i_moved: usize = 0;
        for i_resid in 0..N_RES_MOVED {
            for i_atom in 0..4 {
                v = system.atom_to_nearest_vec3(i_chain, moved_from*4 - 1);
                roto.apply_mut(&mut v);
                proposal.bbx[i_resid] = system.real_to_int( v.x);
                proposal.bby[i_resid] = system.real_to_int(v.y);
                proposal.bbz[i_resid] = system.real_to_int(v.z);
                i_chain += 1;
            }
        }
    }
}