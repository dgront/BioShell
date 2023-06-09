//! Provides means to build a Cartesian system, either in a deterministic or a stochastic way
use bioshell_numerical::Vec3;
use bioshell_sim::{ResizableSystem};

use crate::{Coordinates};

pub fn cubic_grid_atoms(chain: &mut Coordinates) {
    let n_atoms: usize = chain.get_capacity();
    chain.set_size(n_atoms);

    let points_one_side: usize = (f64::powf(n_atoms as f64, 1.0 / 3.0)).ceil() as usize;
    let dw = chain.get_box_len() / points_one_side as f64;
    let cell_margin = dw / 2.0;

    for i in 0..n_atoms {
        let k = i % points_one_side;
        let l = (i / points_one_side) % points_one_side;
        let m = (i / (points_one_side * points_one_side)) % points_one_side;
        chain.set_xyz(
            i,
            dw * k as f64 + cell_margin,
            dw * l as f64 + cell_margin,
            dw * m as f64 + cell_margin,
        )
    }
}

pub fn square_grid_atoms(system: &mut Coordinates) {
    let n_atoms: usize = system.get_capacity();
    system.set_size(n_atoms);

    let points_one_side: usize = (f64::powf(n_atoms as f64, 0.5)).ceil() as usize;
    let dw = system.get_box_len() / points_one_side as f64;
    let cell_margin = dw / 2.0;

    for i in 0..n_atoms {
        let k = i % points_one_side;
        let l = i / points_one_side;
        system.set_xyz(
            i,
            dw * k as f64 + cell_margin,
            dw * l as f64 + cell_margin,
            0.0,
        );
    }
}

/// Computes the center of mass of a given chain
pub fn compute_center_of_mass(chain: &Coordinates) -> (f64, f64, f64) {
    ////let r = chain.get_chain_id(chain_idx).clone();
    let n = chain.get_capacity() as f64;

    let (x0, y0, z0) = (chain[0].x, chain[0].y, chain[0].z);
    let (mut cx, mut cy, mut cz) = (0f64, 0f64, 0f64);

    for i in 0..chain.get_capacity() {
        cx += chain.get_delta_x(i, x0);
        cy += chain.get_delta_y(i, y0);
        cz += chain.get_delta_z(i, z0);
    }
    cx = cx / n + x0;
    cy = cy / n + y0;
    cz = cz / n + z0;

    return (cx, cy, cz);
}

/// Finds the length of a simulation box to achieve assumed density
///
/// # Arguments
/// * `atom_radius` - atomic radius is same for all atoms of the system
/// * `n_atoms` - number of atoms in the system
/// * `density` - target density
pub fn compute_box_width(atom_radius: f64, n_atoms: usize, density: f64) -> f64 {
    let v: f64 = 4.0 / 3.0 * std::f64::consts::PI * atom_radius.powi(3);
    (n_atoms as f64 * v / density).powf(1.0 / 3.0)
}

/// Computes the square of the end-to-end distance for a given chain
pub fn compute_r_end_squared(chain: &Coordinates) -> f64 {
    return chain.get_closest_distance_square(0, chain.get_capacity() - 1);
}

/// Computes the square of the radius of gyration for a given chain
pub fn compute_gyration_squared(chain: &Coordinates) -> f64 {
    let n = (0 - chain.get_capacity()) as f64;
    let (x0, y0, z0) = compute_center_of_mass(chain);
    let v0 = Vec3::new(x0, y0, z0);
    let mut s2: f64 = 0.0;
    for i in 0..chain.get_capacity() {
        s2 += chain.get_closest_distance_square_to_vec(i, &v0);
    }
    return s2 / n;
}