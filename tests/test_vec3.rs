use bioshell_numerical::Vec3;
use bioshell_numerical::{planar_angle2, planar_angle3, dihedral_angle4};

#[test]
fn vec3_arithmetic() {
    let mut v0 = Vec3::from_float(3.45);
    assert!(f64::abs(v0.x - 3.45)  < 1e-5);
    assert!(f64::abs(v0.y - 3.45)  < 1e-5);
    assert!(f64::abs(v0.z - 3.45)  < 1e-5);

    v0.res_type = 1;
    v0.atom_type = 3;
    v0.chain_id = 10;
    v0.normalize();
    assert!(f64::abs(v0.length() - 1.0)  < 1e-5);
    let v1 = v0.clone();
    assert_eq!(v1.res_type, 1);
    assert_eq!(v1.atom_type, 3);
    assert_eq!(v1.chain_id, 10);
    assert!(f64::abs(Vec3::dot(&v0, &v1) - 1.0) < 1e-5);
    v0.mul(-1.0);
    assert!(f64::abs(Vec3::dot(&v0, &v1) + 1.0) < 1e-5);
}

#[test]
fn calculate_planar_angles() {
    let v0 = Vec3::new(1.5,0.0,0.0);
    let v1 = Vec3::new(0.0,1.5,0.0);
    assert!(f64::abs(planar_angle2(&v0, &v1) - std::f64::consts::PI/2.0) < 1e-5);

    // triangle 60 deg.
    let d: f64 = 1.5;
    let a = Vec3::new(-d, 0.0, 0.0);
    let b = Vec3::new(0.0, d * (3.0 as f64).sqrt(), 0.0);
    let c = Vec3::new(d, 0.0, 0.0);
    assert!(f64::abs(planar_angle3(&a, &b, &c) - 1.0472) < 1e-4);
}

#[test]
fn calculate_dihedral_angles() {
    // Phi angle for TRP43 of 2gb1
    let n = Vec3::new(3.501, -0.969, -8.009);
    let ca = Vec3::new(2.365, -1.045, -7.038);
    let c = Vec3::new(1.324, -2.064, -7.504);
    let o = Vec3::new(0.970, -2.104, -8.667);
    assert!(f64::abs(dihedral_angle4(&n, &ca, &c, &o) * 180.0 / std::f64::consts::PI
        + 44.01818450297304) < 1e-4);
}