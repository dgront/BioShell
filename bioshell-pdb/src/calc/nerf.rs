use crate::calc::{Matrix3x3, Vec3};


pub fn place_atom(a: &Vec3, b: &Vec3, c: &Vec3, r: f64, planar: f64, dihedral: f64, output: &mut Vec3) {
    let mut bc = Vec3::sub_s(c, b);
    bc.normalize();
    let mut n = Vec3::sub_s(b,a);
    n.normalize();
    n = Vec3::cross(&n,&bc);
    n.normalize();
    let cross = Vec3::cross(&n, &bc);
    let rot_m = Matrix3x3::from_row_vectors(&bc, &cross, &n);

    let angle = std::f64::consts::PI - planar;

    output.set3(r * angle.cos(),
                r * angle.sin() * dihedral.cos(),
                r * angle.sin() * dihedral.sin());
    println!("befr>>>> {} {} {} {} {} {}", planar.to_degrees(), angle.to_degrees(), dihedral, dihedral.to_degrees(), c.x, output.x);
    rot_m.mul_vec_mut( output);
    println!("{:?}",rot_m);
    println!("aftr>>>> {} {} {} {} {} {}", planar.to_degrees(), angle.to_degrees(), dihedral, dihedral.to_degrees(), c.x, output.x);
    *output += c;
}

/// Rebuilds Cartesian coordinates of a whole chain based on internal coordinates of its atoms.
///
/// i-th atom of a chain is reconstructed based on positions of (i-1), (i-2) and (i-3) atoms,
/// `r[i]` radius, `planar[i]` planar and `dihedral[i]` dihedral; obtained coordinates are writen to `output_chain[i]`.
///
/// # Arguments
///  * `r` - array of bond lengths
///  * `planar` - array of planar angles
///  * `dihedral` - array of dihedral angles in radians
///  * `output_chain` - array of [`Vec3`] objects where the recovered Cartesian coordinates will be stored;
///    **Note**: reconstruction starts from the second atom (i.e. `output_chain[1]`) while
///    `output_chain[0]` is used as the starting point.
///
pub fn place_chain(r: &[f64], planar: &[f64], dihedral: &[f64], output_chain: &mut [Vec3]) {
    if r.len() != planar.len() || r.len() != dihedral.len() || r.len() != output_chain.len() {
        panic!("Inconsistent size of r ({}), alpha ({}), dihedral ({}) and output_chain ({})",
               r.len(), planar.len(), dihedral.len(), output_chain.len())
    }
    output_chain[1].set3(r[1] + output_chain[0].x, output_chain[0].y, output_chain[0].z);
    let angle = std::f64::consts::PI - planar[2];
    output_chain[2].set3(r[2]*angle.cos() + output_chain[1].x,
                         r[2]*angle.sin() + output_chain[1].y, output_chain[1].z );

    for i in 3..r.len() {
        let [a, b, c, d] = &mut output_chain[i-3..i+1] else {
            unreachable!();
        };
        place_atom(a, b, c, r[i], planar[i], dihedral[i], d);
    }
}