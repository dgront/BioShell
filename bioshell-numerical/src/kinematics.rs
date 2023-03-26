use std::f64::consts::PI;
use crate::Vec3;

/// Builds a stub of a molecule by computing coordinates of its first three atoms.
///
/// Cartesian coordinates of the two atoms `j` and `k` of a molecule are computed from its internal degrees of freedom:
/// two bond lengths and a planar angle. Position of the first atom (which is `i`) is used as
/// the starting point for the molecule.
///
/// ```rust
/// use bioshell_numerical::vec3::Vec3;
/// use bioshell_numerical::kinematics::fill_molecule_stub;
/// // This test builds a water molecule from its internal coordinates
/// let mut hoh = vec![Vec3::new(0.0, 0.0, 0.0); 3];
/// fill_molecule_stub(&mut hoh, 0, 1, 2, 0.9485, 0.9485, 104.45_f64.to_radians());
/// assert!((hoh[2].x - 1.18518).abs() < 0.00001);
/// assert!((hoh[2].y - 0.9185).abs() < 0.00001);
/// ```
pub fn fill_molecule_stub(coords: &mut Vec<Vec3>, i: usize, j: usize, k: usize, bond_12: f64, bond_23: f64, angle: f64) {

    coords[j].x = coords[i].x + bond_12;
    coords[j].y = coords[i].y;
    coords[j].z = coords[i].z;
    coords[k].x = coords[j].x + bond_23 * (PI - angle).cos();
    coords[k].y = coords[j].y + bond_23 * (PI - angle).sin();
    coords[k].z = coords[i].z;
}

/// Computes Cartesian coordinates for an atom given its internal coordinates
///
/// Cartesian coordinates of the atoms `l` is evaluated in the reference frame defined by `i`, `j` and `k` atoms
///
/// # Arguments
/// * `i` - index of the first reference atom
/// * `j` - index of the second reference atom
/// * `k` - index of the third reference atom
/// * `l` - index of the reconstructed atom
/// * `bond_length` - distance between atoms `l` and `k`
/// * `planar_angle` - `j`, `l` and `k` planar angle
/// * `dihedral_angle` - `i`, `j`, `l` and `k` dihedral angle
///
/// ```rust
/// use bioshell_numerical::vec3::Vec3;
/// use bioshell_numerical::{dihedral_angle4};
/// use bioshell_numerical::kinematics::{fill_molecule_stub, internal_to_cartesian};
/// // This test builds a methane molecule from its internal coordinates
/// let mut ch4 = vec![Vec3::new(0.0, 0.0, 0.0); 5];
/// let planar_angle = 109.471_f64.to_radians();
/// let dihedral_angle = 120.0_f64.to_radians();
/// fill_molecule_stub(&mut ch4, 0, 1, 2, 1.089, 1.089, planar_angle);
/// internal_to_cartesian(&mut ch4, 0, 2, 1, 3, 1.089, planar_angle, dihedral_angle);
/// internal_to_cartesian(&mut ch4, 0, 2, 1, 4, 1.089, planar_angle, -dihedral_angle);
/// assert!((dihedral_angle4(&ch4[0], &ch4[1], &ch4[2], &ch4[4]).abs() - dihedral_angle).abs() < 0.00001);
/// ```
pub fn internal_to_cartesian(coords: &mut Vec<Vec3>, i: usize, j: usize, k: usize, l: usize,
                             bond_length: f64, planar_angle: f64, dihedral_angle: f64) {

    let sin_planar = (PI - planar_angle).sin();
    let cos_planar = (PI - planar_angle).cos();
    let sin_dih = dihedral_angle.sin();
    let cos_dih = dihedral_angle.cos();
    let bs = bond_length * sin_planar;
    coords[l].x = bond_length * cos_planar;
    coords[l].y = bs * cos_dih;
    coords[l].z = bs * sin_dih;

    // Translate atoms so coords[k] is at the origin.
    let mut a3 = Vec3::sub_two(&coords[i], &coords[k]);
    let a2 = Vec3::sub_two(&coords[j], &coords[k]);

    let d2_xz = a2.x * a2.x + a2.z * a2.z;
    let d2 = (d2_xz + a2.y * a2.y).sqrt();
    let dxz = d2_xz.sqrt();

    let (d2_inverse, dxz_inverse): (f64, f64);
    let (xx1, x2o, y2o, z2o, xz2o): (f64, f64,f64, f64,f64);
    if d2 < 0.0001 {d2_inverse = 1.0 / 0.001;}
    else {d2_inverse = 1.0 / d2;}

    if dxz < 0.001 {
        xx1 = a3.x;
        x2o = 1.0;
        z2o = 0.0;
    } else {
        dxz_inverse = 1.0 / dxz;
        x2o = a2.x * dxz_inverse;
        z2o = a2.z * dxz_inverse;
        xx1 = a3.x * x2o + a3.z * z2o;
        a3.z = a3.z * x2o - a3.x * z2o;
    }

    xz2o = dxz * d2_inverse;
    y2o = a2.y * d2_inverse;
    a3.x = (-xx1 * xz2o) - a3.y * y2o;
    a3.y = xx1 * y2o - a3.y * xz2o;

    let dyz = (a3.y * a3.y + a3.z * a3.z).sqrt();
    let dyz_inverse = 1.0 / dyz;
    let y1o = a3.y * dyz_inverse;
    let z1o = a3.z * dyz_inverse;
    let yy4 = y1o * coords[l].y - z1o * coords[l].z;
    let zz4 = y1o * coords[l].z + z1o * coords[l].y;

    let xx4 = y2o * yy4 - xz2o * coords[l].x;
    coords[l].y = (-xz2o * yy4) - y2o * coords[l].x;
    coords[l].x = x2o * xx4 - z2o * zz4;
    coords[l].z = z2o * xx4 + x2o * zz4;

    coords[l].x += coords[k].x;
    coords[l].y += coords[k].y;
    coords[l].z += coords[k].z;
}