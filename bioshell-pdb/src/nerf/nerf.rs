use crate::calc::{Matrix3x3, Vec3};

pub fn create_stub(r: &[f64], planar: &[f64], output_chain: &mut [Vec3]) {

    if r.len() != 3|| planar.len()  != 3 || output_chain.len() != 3 {
        panic!("build_stub() function requires 3-element slices for distances, angles and positions")
    }
    output_chain[1].set3(r[1] + output_chain[0].x, output_chain[0].y, output_chain[0].z);
    let angle = std::f64::consts::PI - planar[2];
    output_chain[2].set3(r[2]*angle.cos() + output_chain[1].x,
                         r[2]*angle.sin() + output_chain[1].y, output_chain[1].z );
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


/// Rebuilds Cartesian coordinates of a linear chain given internal coordinates of its atoms.
///
/// This simplified function reconstructs i-th atom of a chain  based on positions of (i-1), (i-2) and (i-3) atoms,
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

    create_stub(&r[0..3], &planar[0..3], &mut output_chain[0..3]);

    for i in 3..r.len() {
        let [a, b, c, d] = &mut output_chain[i-3..i+1] else {
            unreachable!();
        };
        restore_atom(a, b, c, r[i], planar[i], dihedral[i], d);
    }
}

pub fn restore_branched_chain(r: &[f64], planar: &[f64], dihedral: &[f64], topology: &[[usize;4]],
                              output_chain: &mut [Vec3]) {

    if r.len() != planar.len() || r.len() != dihedral.len() || r.len() != output_chain.len() {
        panic!("Inconsistent size of r ({}), alpha ({}), dihedral ({}) and output_chain ({})",
               r.len(), planar.len(), dihedral.len(), output_chain.len())
    }

    // --- rebuild the stub
    create_stub(&r[0..3], &planar[0..3], &mut output_chain[0..3]);

    let mut v = Vec3::default();
    for iatom in 3..output_chain.len() {
        let [i, j, k, l] = topology[iatom];
        restore_atom(&output_chain[i], &output_chain[j], &output_chain[k],
                     r[iatom], planar[iatom], dihedral[iatom], &mut v);
        output_chain[l].set(&v);
    }

}