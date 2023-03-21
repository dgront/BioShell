use bioshell_numerical::Vec3;

use crate::Coordinates;

/// Computes the center of mass of a given chain
pub fn cm(chain: &Coordinates, chain_idx: usize) -> (f64, f64, f64) {
    let r = chain.chain_range(chain_idx).clone();
    let n = (r.end - r.start) as f64;

    let (x0, y0, z0) = (chain[r.start].x, chain[r.start].y, chain[r.start].z);
    let (mut cx, mut cy, mut cz) = (0f64, 0f64, 0f64);
    for ai in r {
        cx += chain.delta_x(ai, x0);
        cy += chain.delta_y(ai, y0);
        cz += chain.delta_z(ai, z0);
    }
    cx = cx / n + x0;
    cy = cy / n + y0;
    cz = cz / n + z0;
    return (cx, cy, cz);
}

/// Computes the square of the end-to-end distance for a given chain
pub fn r_end_squared(chain: &Coordinates, chain_idx: usize) -> f64 {
    let chain_range = chain.chain_range(chain_idx);

    return chain.closest_distance_square(chain_range.start, chain_range.end - 1);
}

/// Computes the square of the radius of gyration for a given chain
pub fn gyration_squared(chain: &Coordinates, chain_idx: usize) -> f64 {
    let r = chain.chain_range(chain_idx).clone();
    let n = (r.end - r.start) as f64;
    let (x0, y0, z0) = cm(chain, chain_idx);
    let v0 = Vec3::new(x0, y0, z0);
    let mut s2: f64 = 0.0;
    for ai in r {
        s2 += chain.closest_distance_square_to_vec(ai, &v0);
    }

    return s2 / n;
}
