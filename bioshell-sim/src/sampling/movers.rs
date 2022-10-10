use std::ops::Range;
use rand::Rng;

use bioshell_ff::{System};


/// performs a volume change
#[allow(non_snake_case)]
pub fn change_volume(system: &mut System, max_step:f64) -> Range<usize> {
    let mut rng = rand::thread_rng();

    let v0 = system.volume();
    let lnV0 = v0.ln();
    let lnV = lnV0 + rng.gen_range(-max_step..max_step);
    let new_len = lnV.exp().powf(0.333333333);
    system.set_box_len(new_len);

    0..system.size()
}

pub fn single_atom_move(future: &mut System, max_step:f64) -> Range<usize> {
    let mut rng = rand::thread_rng();
    let i_moved = rng.gen_range(0..future.size());
    future.add(i_moved,rng.gen_range(-max_step..max_step),
               rng.gen_range(-max_step..max_step),rng.gen_range(-max_step..max_step));

    i_moved..i_moved
}

pub fn perturb_chain_fragment(chains: &mut System, max_step:f64) -> Range<usize> {

    const N: usize = 3;
    const F: f64 = 2.0 / (1.0 + N as f64);

    let mut rng = rand::thread_rng();
    let mut moved_from = rng.gen_range(0..chains.size()-N);
    let mut moved_to = moved_from + N - 1;
    while chains.coordinates()[moved_from].chain_id != chains.coordinates()[moved_to].chain_id {
        moved_from = rng.gen_range(0..chains.size() - N);
        moved_to = moved_from + N;
    }

    let dx: f64 = rng.gen_range(-max_step..max_step) * F;
    let dy: f64 = rng.gen_range(-max_step..max_step) * F;
    let dz: f64 = rng.gen_range(-max_step..max_step) * F;

    for i in 0..N/2 {
        let fi: f64 = (i + 1) as f64;
        chains.add(moved_from + i, dx * fi, dy * fi, dz * fi);
        chains.add(moved_to - i, dx * fi, dy * fi, dz * fi);
    }

    if N % 2 == 1 {
        let mid = (moved_from + moved_to) / 2;
        chains.add(mid,dx / F, dy / F, dz / F);
    }

    moved_from..moved_to
}


