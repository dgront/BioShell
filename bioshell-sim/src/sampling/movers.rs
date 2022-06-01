use std::ops::Range;
use rand::Rng;

use bioshell_ff::Coordinates;


pub fn single_atom_move(future: &mut Coordinates, max_step:f32) -> Range<usize> {
    let mut rng = rand::thread_rng();
    let i_moved = rng.gen_range(0..future.size());
    future[i_moved].x += rng.gen_range(-max_step..max_step);
    future[i_moved].z += rng.gen_range(-max_step..max_step);
    future[i_moved].y += rng.gen_range(-max_step..max_step);
    i_moved..i_moved
}

pub fn perturb_chain_fragment() {}


