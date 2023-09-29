use rand::{Rng, thread_rng};
use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::{MoveProposal, Mover, SurpassAlphaSystem};
use crate::tail_move::MovedTermini::{CTerminal, NTerminal};

pub struct TailMove {
    max_angle: f64,
    max_range_allowed: f64
}

pub enum MovedTermini {
    NTerminal,
    CTerminal
}

impl Mover<1> for TailMove {
    fn propose(&self, system: &SurpassAlphaSystem, proposal: &mut MoveProposal<1>) {

        let mut rng = thread_rng();
        // --- pick a chain randomly
        let i_chain = rng.gen_range(0..system.count_chains());
        // --- pick either end
        let which_tail = match rng.gen::<bool>() {
            true => {NTerminal}
            false => {CTerminal}
        };
        let angle = rng.gen_range(-self.max_angle..self.max_angle);

        self.compute_move(system, i_chain, which_tail, angle, proposal);
    }
}

impl TailMove {

    pub fn new(max_angle: f64, max_range_allowed: f64) -> TailMove {
        TailMove{ max_angle, max_range_allowed }
    }

    pub fn compute_move(&self, system: &SurpassAlphaSystem, which_chain: usize, which_tail: MovedTermini,
                        angle: f64, proposal: &mut MoveProposal<1>) {

        let r = system.chain_atoms(which_chain);
        let (axis_from, axis_to, mut moved, i_moved, i_referenced) = match which_tail {
            CTerminal => {
                (system.ca_to_nearest_vec3(r.end-3, r.end-2),
                 system.ca_to_vec3(r.end-2),
                 system.ca_to_nearest_vec3(r.end-1, r.end-2),
                 r.end-1, r.end - 2)
            }
            NTerminal => {
                (system.ca_to_vec3(r.start + 1),
                 system.ca_to_nearest_vec3(r.start+2, r.start+1),
                 system.ca_to_nearest_vec3(r.start, r.start+1),
                 r.start, r.start + 1)
            }
        };
        // --- prepare rototranslation
        let roto = Rototranslation::around_axis(&axis_from, &axis_to, angle);

        roto.apply_mut(&mut moved);
        proposal.cax[0] = system.real_to_int(moved.x);
        proposal.cay[0] = system.real_to_int(moved.y);
        proposal.caz[0] = system.real_to_int(moved.z);
        proposal.first_moved_pos = i_moved;
    }
}