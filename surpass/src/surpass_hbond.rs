use bioshell_pdb::calc::{Rototranslation, Vec3};

/// Calculates x, r, alpha coordinates of a CA atom from its Cartesian coordinates
///
/// This function computes internal coordinates of the ``ca_j`` atom in the Local coordinate System
/// formed by the three atoms: ``ca_prev``, ``ca_i`` and ``ca_next``. Subsequently it computes
/// the `x``, ``r`` and ``alpha`` corrected coordinates of that point
fn compute_x_r_a(ca_prev: &Vec3, ca_i: &Vec3, ca_next: &Vec3, ca_j: &Vec3,
                 correct_y: f64, correct_z: f64) -> (f64, f64, f64) {

    let rt = Rototranslation::by_three_atoms(&ca_prev, &ca_i, &ca_next);
    let mut ca_j_local = rt.apply(&ca_j);

    let x = ca_j_local.x;                                        // --- X coordinate
    ca_j_local.y -= correct_y;                                        // --- correct for the center of Y-Z circle
    ca_j_local.z -= correct_z;
    // if (ca_j.z > 4.0) return false;                          // --- "stacking" strands is not allowed!
    let r = (ca_j_local.y * ca_j_local.y + ca_j_local.z * ca_j_local.z).sqrt();    // --- polar coordinates from Y-Z
    let mut a = ca_j_local.z.atan2(ca_j_local.y);                // --- angle should be [-pi,pi]
    if a < 0.0 { a += 2 * std::f64::consts::PI; }               // --- take care the angle is in the range [0,2pi]
    if a < std::f64::consts::PI / 2.0 { a += 2 * std::f64::consts::PI; }

    return (x, r, a);
}