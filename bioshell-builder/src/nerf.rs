//! Implementation of the Natural Reference Frame (NeRF) algorithm that is used to reconstruct Cartesian positions from internal coordinates
//!
use bioshell_pdb::calc::{Matrix3x3, Vec3};

/// Creates a three-atom stub where a fourth atom may be placed into.
pub fn create_stub(a: &Vec3, r_ab: f64, r_bc: f64, a_abc: f64, b: &mut Vec3, c: &mut Vec3) {

    b.set3(a.x+r_ab,a.y, a.z);
    let angle = std::f64::consts::PI - a_abc;
    c.set3(r_bc*angle.cos() + b.x, r_bc*angle.sin() + b.y, b.z );
}

pub fn restore_atom(a: &Vec3, b: &Vec3, c: &Vec3, r: f64, planar: f64, dihedral: f64, output: &mut Vec3) {
    let mut bc = Vec3::sub_s(c, b);
    bc.normalize();
    let mut n = Vec3::sub_s(b,a);
    n.normalize();
    n = Vec3::cross(&n,&bc);
    n.normalize();
    let cross = Vec3::cross(&n, &bc);
    let rot_m = Matrix3x3::from_column_vectors(&bc, &cross, &n);

    let angle = std::f64::consts::PI - planar;

    output.set3(r * angle.cos(),
                r * angle.sin() * dihedral.cos(),
                r * angle.sin() * dihedral.sin());
    rot_m.mul_vec_mut( output);
    *output += c;
}


/// Calculates Cartesian coordinates of a linear chain given internal coordinates of its atoms.
///
/// This simplified function reconstructs every i-th atom of a chain  based on positions of (i-1), (i-2) and (i-3) atoms,
/// `r[i]` radius, `planar[i]` planar and `dihedral[i]` dihedral; resulting coordinates are writen to `output_chain[i]`.
///
/// # Arguments
///  * `r` - array of bond lengths
///  * `planar` - array of planar angles
///  * `dihedral` - array of dihedral angles in radians
///  * `output_chain` - array of [`Vec3`] objects where the recovered Cartesian coordinates will be stored;
///    **Note**: reconstruction starts from the second atom (i.e. `output_chain[1]`) while
///    `output_chain[0]` is used as the starting point.
///
pub fn restore_linear_chain(r: &[f64], planar: &[f64], dihedral: &[f64], output_chain: &mut [Vec3]) {

    if r.len() != planar.len() || r.len() != dihedral.len() || r.len() != output_chain.len() {
        panic!("Inconsistent size of r ({}), alpha ({}), dihedral ({}) and output_chain ({})",
               r.len(), planar.len(), dihedral.len(), output_chain.len())
    }

    let [a, b, c] = &mut output_chain[0..3] else { unreachable!(); };
    create_stub(&a, r[1], r[2], planar[2], b,  c);

    for i in 3..r.len() {
        let [a, b, c, d] = &mut output_chain[i-3..i+1] else { unreachable!(); };
        restore_atom(a, b, c, r[i], planar[i], dihedral[i], d);
    }
}

/// Calculates Cartesian coordinates of a molecule from internal coordinates of its atoms.
///
/// ```
/// # use bioshell_builder::nerf::restore_branched_chain;
/// use bioshell_pdb::calc::Vec3;
/// let r_CH: f64 = 1.05;
/// let a_HCH: f64 = 109.471_f64.to_radians();
/// let r = vec![0.0, r_CH, r_CH, r_CH, r_CH];
/// let a = vec![0.0, 0.0, a_HCH, a_HCH, a_HCH];
/// let t = vec![0.0, 0.0, 0.0, 120.0_f64.to_radians(), 240.0_f64.to_radians()];
/// let topo = vec![[0, 0, 0, 0], [0, 1, 0, 0], [1, 0, 2, 0], [1, 2, 0, 3], [1, 2, 0, 4]];
/// let mut methane = vec![Vec3::default(); 5];
/// restore_branched_chain(&r, &a, &t, &topo, &mut methane);
/// ```
/// # Arguments
///  * `r` - array of bond lengths
///  * `planar` - array of planar angles
///  * `dihedral` - array of dihedral angles in radians
///  * `topology` - array that provides indexes of atoms
///  * `output_chain` - array of [`Vec3`] objects where the recovered Cartesian coordinates will be stored;
///    **Note**: reconstruction starts from the second atom (i.e. `output_chain[1]`) while
///    `output_chain[0]` is used as the starting point.
///
pub fn restore_branched_chain(r: &[f64], planar: &[f64], dihedral: &[f64], topology: &[[usize;4]],
                              output_chain: &mut [Vec3]) {

    if r.len() != planar.len() || r.len() != dihedral.len() || r.len() != output_chain.len() {
        panic!("Inconsistent size of r ({}), alpha ({}), dihedral ({}) and output_chain ({})",
               r.len(), planar.len(), dihedral.len(), output_chain.len())
    }

    // --- rebuild the stub
    let mut b = Vec3::default();
    let mut v = Vec3::default();
    let [k, l, m] = topology[2][1..4] else { unreachable!(); };
    create_stub(&output_chain[k], r[1], r[2], planar[2], &mut b, &mut v);
    output_chain[l].set(&b);
    output_chain[m].set(&v);

    for iatom in 3..output_chain.len() {
        let [i, j, k, l] = topology[iatom];
        restore_atom(&output_chain[i], &output_chain[j], &output_chain[k],
                     r[iatom], planar[iatom], dihedral[iatom], &mut v);
        output_chain[l].set(&v);
    }
}

/// Calculates Cartesian coordinates of a molecule from internal coordinates of its atoms.
///
/// This function allows restoring atoms in a user-defined order
/// # Arguments
///  * `r` - array of bond lengths
///  * `planar` - array of planar angles
///  * `dihedral` - array of dihedral angles in radians
///  * `topology` - array that provides indexes of atoms
///  * `output_chain` - array of [`Vec3`] objects where the recovered Cartesian coordinates will be stored;
///    **Note**: reconstruction starts from the second atom (i.e. `output_chain[1]`) while
///    `output_chain[0]` is used as the starting point.
///  * `restoring_order` - order in which Cartesian coordinates evaluated
///
pub fn restore_branched_chain_in_order(r: &[f64], planar: &[f64], dihedral: &[f64], topology: &[[usize;4]],
                              restoring_order: &Vec<usize>, output_chain: &mut [Vec3]) {

    if r.len() != planar.len() || r.len() != dihedral.len() || r.len() != output_chain.len() {
        panic!("Inconsistent size of r ({}), alpha ({}), dihedral ({}) and output_chain ({})",
               r.len(), planar.len(), dihedral.len(), output_chain.len())
    }

    // --- rebuild the stub
    let mut b = Vec3::default();
    let mut v = Vec3::default();
    let [k, l, m] = topology[2][1..4] else { unreachable!(); };
    create_stub(&output_chain[k], r[1], r[2], planar[2], &mut b, &mut v);
    output_chain[l].set(&b);
    output_chain[m].set(&v);

    for index in 3..output_chain.len() {
        let iatom = restoring_order[index];
        let [i, j, k, l] = topology[iatom];
        restore_atom(&output_chain[i], &output_chain[j], &output_chain[k],
                     r[iatom], planar[iatom], dihedral[iatom], &mut v);
        output_chain[l].set(&v);
    }
}