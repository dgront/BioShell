use bioshell_numerical::Vec3;
use bioshell_numerical::kinematics::{fill_molecule_stub, internal_to_cartesian};
use bioshell_numerical::{dihedral_angle4, planar_angle3};

#[test]
fn build_water() {
    let mut hoh = vec![Vec3::new(0.0, 0.0, 0.0); 3];
    fill_molecule_stub(&mut hoh, 0, 1, 2, 0.9485, 0.9485, 104.45_f64.to_radians());
    assert!((hoh[1].x - 0.9485).abs() < 0.00001);
    assert!((hoh[1].y).abs() < 0.00001);
    assert!((hoh[2].x - 1.18518).abs() < 0.00001);
    assert!((hoh[2].y - 0.9185).abs() < 0.00001);
    let mut v = Vec3::new(0.0, 1.0, 2.0);
    v.res_type = 1;
    v.atom_type = 6;
    v.chain_id = 128;
    println!("{:?}", v);
}

#[test]
fn build_methane() {
    let mut ch4 = vec![Vec3::new(0.0, 0.0, 0.0); 5];
    let planar_angle = 109.471_f64.to_radians();
    let dihedral_angle = 120.0_f64.to_radians();

    fill_molecule_stub(&mut ch4, 0, 1, 2, 1.089, 1.089, planar_angle);
    internal_to_cartesian(&mut ch4, 0, 2, 1, 3, 1.089, planar_angle, dihedral_angle);
    internal_to_cartesian(&mut ch4, 0, 2, 1, 4, 1.089, planar_angle, -dihedral_angle);

    let planars = vec![(0, 1, 2), (0, 1, 3), (0, 1, 4), (2, 1, 3), (2, 1, 4)];
    for (i, j, k) in planars {
        assert!((planar_angle3(&ch4[i], &ch4[j], &ch4[k]) - planar_angle).abs() < 0.00001);
    }
    let dihedrals = vec![(0, 1, 2, 3), (0, 1, 3, 4), (2, 1, 3, 4)];
    for (i, j, k, l) in dihedrals {
        assert!((dihedral_angle4(&ch4[i], &ch4[j], &ch4[k], &ch4[l]).abs() - dihedral_angle).abs() < 0.00001);
    }
}