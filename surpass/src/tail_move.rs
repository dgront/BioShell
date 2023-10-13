use std::fmt;
use log::debug;
use rand::{Rng};
use bioshell_pdb::calc::{Rototranslation, Vec3};
use crate::{MoveProposal, Mover, SurpassAlphaSystem};
use crate::tail_move::MovedTermini::{CTerminal, NTerminal};

pub enum MovedTermini {
    NTerminal,
    CTerminal
}

impl fmt::Display for MovedTermini {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            MovedTermini::NTerminal => write!(f, "NTerminal"),
            MovedTermini::CTerminal => write!(f, "CTerminal"),
        }
    }
}

pub struct TailMove<const N_MOVED: usize> {
    max_angle: f64,
    max_range_allowed: f64
}

impl<const N_MOVED: usize> Mover<N_MOVED> for TailMove<N_MOVED> {
    fn propose<R: Rng>(&self, system: &SurpassAlphaSystem, rnd_gen: &mut R, proposal: &mut MoveProposal<N_MOVED>) {

        // --- pick a chain randomly
        let i_chain = rnd_gen.gen_range(0..system.count_chains());
        // --- pick either end
        let which_tail = match rnd_gen.gen::<bool>() {
            true => {NTerminal}
            false => {CTerminal}
        };
        let angle = rnd_gen.gen_range(-self.max_angle..self.max_angle);

        self.compute_move(system, i_chain, which_tail, angle, proposal);
    }
}

impl<const N_MOVED: usize> TailMove<N_MOVED> {

    pub fn new(max_angle: f64, max_range_allowed: f64) -> TailMove<N_MOVED> {
        TailMove{ max_angle, max_range_allowed }
    }

    pub fn compute_move(&self, system: &SurpassAlphaSystem, which_chain: usize, which_tail: MovedTermini,
                        angle: f64, proposal: &mut MoveProposal<N_MOVED>) {

        let r = system.chain_atoms(which_chain);
        let (axis_from, axis_to, mut moved, i_moved, i_referenced) = match which_tail {
            CTerminal => {
                (system.ca_to_nearest_vec3(r.end-N_MOVED-2, r.end-N_MOVED-1),
                 system.ca_to_vec3(r.end-2),
                 system.ca_to_nearest_vec3(r.end-N_MOVED, r.end-N_MOVED-1),
                 r.end-N_MOVED, r.end - N_MOVED -1)
            }
            NTerminal => {
                (system.ca_to_vec3(r.start + N_MOVED),
                 system.ca_to_nearest_vec3(r.start+N_MOVED+1, r.start+N_MOVED),
                 system.ca_to_nearest_vec3(r.start, r.start+N_MOVED),
                 r.start, r.start + N_MOVED)
            }
        };
        debug!("{} tail move of {}", which_tail, i_moved);

        // --- prepare rototranslation
        let roto = Rototranslation::around_axis(&axis_from, &axis_to, angle);

        roto.apply_mut(&mut moved);
        proposal.cax[0] = system.real_to_int(moved.x);
        proposal.cay[0] = system.real_to_int(moved.y);
        proposal.caz[0] = system.real_to_int(moved.z);
        proposal.first_moved_pos = i_moved;
    }
}