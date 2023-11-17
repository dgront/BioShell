use std::fmt;
use log::trace;
use rand::{Rng};
use bioshell_pdb::calc::{Rototranslation};
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

        let r = system.chain_residues(which_chain);
        // --- axis_from, axis_to : vectors defining the rotation axis
        // --- i_moved_from, i_moved_to : range of the atoms being moved: i_moved_from..i_moved_to
        // --- i_referenced index of the atom used as the PBC reference, i.e. all other atoms are moved to the image which is the closest to i_referenced
        let (axis_from, axis_to, i_moved_from, i_moved_to, i_referenced) = match which_tail {
            CTerminal => {
                (system.atom_to_nearest_vec3(r.end-N_MOVED-2, r.end-N_MOVED-1),
                 system.atom_to_vec3(r.end-N_MOVED-1),
                 r.end - N_MOVED, r.end, r.end - N_MOVED -1)
            }
            NTerminal => {
                (system.atom_to_vec3(r.start + N_MOVED),
                 system.atom_to_nearest_vec3(r.start+N_MOVED+1, r.start+N_MOVED),
                 r.start, r.start + N_MOVED, r.start + N_MOVED)
            }
        };
        trace!("{} tail move of {}..{}", which_tail, i_moved_from, i_moved_to);

        // --- prepare rototranslation
        let roto = Rototranslation::around_axis(&axis_from, &axis_to, angle);

        let mut i_pos = 0;
        for i_moved in i_moved_from..i_moved_to {
            let mut moved = system.atom_to_nearest_vec3(i_moved, i_referenced);
            roto.apply_mut(&mut moved);
            proposal.cax[i_pos] = system.real_to_int(moved.x);
            proposal.cay[i_pos] = system.real_to_int(moved.y);
            proposal.caz[i_pos] = system.real_to_int(moved.z);
            i_pos += 1;
        }
        proposal.first_moved_pos = i_moved_from;
    }
}