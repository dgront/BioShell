macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}

#[cfg(test)]
mod rototranslation_test {

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

        println!("{:?}", rot.rotation_matrix());
        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        for _ in 0..4 {
            rot.apply_mut(&mut another_vec);
            println!("{}", another_vec);
        }
        assert_delta!(another_vec.x, 0.0, 0.000001);
        assert_delta!(another_vec.y, 1.0, 0.000001);
        assert_delta!(another_vec.z, 0.0, 0.000001);
    }

    #[test]
    fn rotate_cube_around_axis_001() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(0.0, 0.0, 1.0);
        let angle = std::f64::consts::PI / 2.0;
        let rot = Rototranslation::around_axis(&corner_1, &corner_2, angle);
        println!("{}", corner_1);
        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        for _ in 0..4 {
            rot.apply_mut(&mut another_vec);
            println!("{}", another_vec);
        }
    }
}
