

#[cfg(test)]
mod rototranslation_test {
    use bioshell_pdb::{assert_delta, assert_vec3_eq};
    use bioshell_pdb::calc::{Vec3, Rototranslation};

    #[test]
    fn rototranslation_apply_mut() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = 180.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(&_start_vector, &_end_vector, _theta_rad);

        let mut _candidate_vector: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        _roto.apply_mut(&mut _candidate_vector);

        assert_delta!(_candidate_vector.x, 3.0, 0.00001);
        assert_delta!(_candidate_vector.y, 2.0, 0.00001);
        assert_delta!(_candidate_vector.z, 1.0, 0.00001);
    }

    #[test]
    fn rototranslation_apply_mut2() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = - 180.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(&_end_vector,&_start_vector,  _theta_rad);

        let mut _candidate_vector: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        _roto.apply_mut(&mut _candidate_vector);

        assert_delta!(_candidate_vector.x, 3.0, 0.00001);
        assert_delta!(_candidate_vector.y, 2.0, 0.00001);
        assert_delta!(_candidate_vector.z, 1.0, 0.00001);
    }

    #[test]
    fn rototranslation_apply_inverse_mut() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = 180.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(&_start_vector, &_end_vector, _theta_rad);

        let mut _candidate_vector: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        _roto.apply_inverse_mut(&mut _candidate_vector);

        assert_delta!(_candidate_vector.x, 3.0, 0.00001);
        assert_delta!(_candidate_vector.y, 2.0, 0.00001);
        assert_delta!(_candidate_vector.z, 1.0, 0.00001);
    }

    #[test]
    fn rotate_cube_around_axis_111() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(1.0, 1.0, 1.0);
        let angle = 2.0 * std::f64::consts::PI / 3.0;
        let rot = Rototranslation::around_axis(&corner_1, &corner_2, angle);

        let expected = [0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 1.00, 0.00];
        for i in 0..9 {
            assert_delta!(rot.rotation_matrix()[i], expected[i], 0.000001);
        }

        // println!("{:?}", rot.rotation_matrix());
        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        for _ in 0..4 {
            rot.apply_mut(&mut another_vec);
            // println!("{}", another_vec);
        }
        assert_delta!(another_vec.x, 0.0, 0.000001, "bad X coordinate after 111 rotation");
        assert_delta!(another_vec.y, 1.0, 0.000001, "bad Y coordinate after 111 rotation");
        assert_delta!(another_vec.z, 0.0, 0.000001, "bad Z coordinate after 111 rotation");
    }

    #[test]
    fn rotate_cube_around_axis_001() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(0.0, 0.0, 1.0);
        let angle = std::f64::consts::PI / 2.0;
        let rot = Rototranslation::around_axis(&corner_1, &corner_2, angle);
        // println!("{}", corner_1);
        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        for _ in 0..4 {
            rot.apply_mut(&mut another_vec);
            // println!("{}", another_vec);
        }
    }

    #[test]
    fn point_in_lcs() {
        let prev_p = Vec3::new(-3.0, -1.0, 0.0);
        let the_p = Vec3::new(0.0, 0.0, 0.0);
        let next_p = Vec3::new(3.0, -1.0, 0.0);
        let rot = Rototranslation::by_three_atoms(&prev_p, &the_p, &next_p);
        let result = format!("{}", rot);
        let expected = r#"[  1.000   0.000   0.000 ]     | vx -  0.000 |
[  0.000   0.000   1.000 ]  *  | vy -  0.000 |
[  0.000  -1.000   0.000 ]     | vz -  0.000 |
"#;
        assert_eq!(result, expected);
        // println!("{}", rot);
        let mut p = Vec3::new(0.0, 0.0, -3.0);
        rot.apply_mut(&mut p);
        let expected = Vec3::new(0.0, -3.0, 0.0);
        assert_vec3_eq!(p, expected, 0.000001, "incorrect vector");
    }
}
