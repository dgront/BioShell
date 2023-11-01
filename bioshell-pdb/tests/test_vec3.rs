
#[cfg(test)]
mod test_vec3 {
    use bioshell_pdb::assert_delta;
    use bioshell_pdb::calc::{Vec3, dihedral_angle4};

    #[test]
    fn test_dihedrals() {

        // some atom coordinates for testing
        let p0 = Vec3::from_array(&[24.969, 13.428, 30.692]); // N
        let p1 = Vec3::from_array(&[24.044, 12.661, 29.808]); // CA
        let p2 = Vec3::from_array(&[22.785, 13.482, 29.543]); // C
        let p3 = Vec3::from_array(&[21.951, 13.670, 30.431]); // O
        let p4 = Vec3::from_array(&[23.672, 11.328, 30.466]); // CB
        let p5 = Vec3::from_array(&[22.881, 10.326, 29.620]); // CG
        let p6 = Vec3::from_array(&[23.691,  9.935, 28.389]); // CD1
        let p7 = Vec3::from_array(&[22.557,  9.096, 30.459]); // CD2

        assert_delta!(dihedral_angle4(&p0, &p1, &p2, &p3).to_degrees(), -71.21515, 1E-4);
        assert_delta!(dihedral_angle4(&p0, &p1, &p4, &p5).to_degrees(), -171.94319, 1E-4);
        assert_delta!(dihedral_angle4(&p1, &p4, &p5, &p6).to_degrees(), 60.82226, 1E-4);
        assert_delta!(dihedral_angle4(&p1, &p4, &p5, &p7).to_degrees(), -177.63641, 1E-4);

        // Check if the sign is correct, i.e. if it follows the UPAC rule
        let a = Vec3::from_array(&[-0.5, 0.0, 1.0]);
        let b = Vec3::from_float(0.0);
        let c = Vec3::from_array(&[0.0, 1.0, 0.0]);
        let d = Vec3::from_array(&[0.0, 1.0, 1.0]);
        assert!(dihedral_angle4(&a, &b, &c, &d) > 0.0);
    }
}