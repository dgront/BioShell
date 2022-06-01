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

pub fn perturb_chain_fragment(chains: &mut Coordinates, max_step:f32) -> Range<usize> {

    const N: usize = 3;
    const F: f32 = 2.0 / (1.0 + N as f32);

    let mut rng = rand::thread_rng();
    let mut moved_from = rng.gen_range(0..chains.size()-N);
    let mut moved_to = moved_from + N - 1;
    while chains[moved_from].chain_id != chains[moved_to].chain_id {
        moved_from = rng.gen_range(0..chains.size() - N);
        moved_to = moved_from + N;
    }

    let dx: f32 = rng.gen_range(-max_step..max_step) * F;
    let dy: f32 = rng.gen_range(-max_step..max_step) * F;
    let dz: f32 = rng.gen_range(-max_step..max_step) * F;

    for i in 0..N/2 {
        let fi: f32 = (i + 1) as f32;
        chains[moved_from + i].add3(dx * fi, dy * fi, dz * fi);
        chains[moved_to - i].add3(dx * fi, dy * fi, dz * fi);
    }

    if N % 2 == 1 {
        let mid = (moved_from + moved_to) / 2;
        chains[mid].add3(dx / F, dy / F, dz / F);
    }

    moved_from..moved_to
}


