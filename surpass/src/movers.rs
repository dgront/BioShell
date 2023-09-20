use std::ops::Range;
use rand::{Rng, thread_rng};
use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::SurpassAlphaSystem;

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


pub struct HingeMove<const HINGE_MOVE_SIZE: usize> {
    max_angle: f64,
    max_range_allowed: f64
}

impl<const HINGE_MOVE_SIZE: usize> HingeMove<HINGE_MOVE_SIZE> {
    pub fn new(max_angle: f64, max_range_allowed: f64) -> HingeMove<HINGE_MOVE_SIZE> {
        HingeMove{ max_angle, max_range_allowed }
    }

    pub fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<HINGE_MOVE_SIZE>) {

        let mut rng = thread_rng();
        // --- pick end points randomly
        let mut moved_from = 0;
        let mut moved_to = HINGE_MOVE_SIZE + 1;
        loop {
            moved_from = rng.gen_range(0..system.count_atoms() - HINGE_MOVE_SIZE - 1);
            moved_to = moved_from + HINGE_MOVE_SIZE - 1;
            if system.chain(moved_from) == system.chain(moved_to) { break };
        }
        let angle = rng.gen_range(-self.max_angle..self.max_angle);

        self.compute_move(system, moved_from, angle, proposal);
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
    pub fn compute_move(&self, system: &SurpassAlphaSystem, moved_from: usize, angle: f64, proposal: &mut MoveProposal<HINGE_MOVE_SIZE>) {

        // -------------------------
        // Moving 8 atoms indexed from 1 to 8 right-exclusive (i.e. moved_from = 1) assumes:
        //  - whole window is 10 atoms long
        //  - pivot positions are 0 and 9
        //  - rotated range is 1..9
        // -------------------------
        // --- prepare rototranslation
        let start_vector: Vec3 = system.ca_to_vec3(moved_from - 1);
        let end_vector: Vec3 = system.ca_to_vec3(moved_from + HINGE_MOVE_SIZE);
        let roto = Rototranslation::around_axis(&start_vector, &end_vector, angle);
        // --- rotate atoms around the axis and copy outcome to a proposal
        proposal.first_moved_pos = moved_from;
        let mut v= Vec3::default();
        let mut i_chain = moved_from;
        for i_moved in 0..HINGE_MOVE_SIZE {
            v.set3(system.int_to_real(system.cax[i_chain]),
                   system.int_to_real(system.cay[i_chain]),
                   system.int_to_real(system.caz[i_chain]));
            roto.apply_mut(&mut v);
            proposal.moved_cax[i_moved] = system.real_to_int( v.x);
            proposal.moved_cay[i_moved] = system.real_to_int(v.y);
            proposal.moved_caz[i_moved] = system.real_to_int(v.z);
            i_chain += 1;
        }
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