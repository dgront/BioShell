#[cfg(test)]
mod rototranslation_test {
    use rand_distr::num_traits::float::FloatCore;
    use bioshell_numerical::{Rototranslation, Vec3, Matrix3x3};

    fn truncate(float: f64)->f64
    {
        return f64::trunc(float  * 100.0) / 100.0; // or f32::trunc
    }

    #[test]
    fn rototranslation_apply_mut() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = 180.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(_origin_vector, _start_vector, _end_vector, _theta_rad);

        let _candidate_vector: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        let _rotated_vec = _roto.apply_mut(&_candidate_vector);

        assert_eq!(truncate(_rotated_vec.x), 3.0);
        assert_eq!(truncate(_rotated_vec.y), 2.0);
        assert_eq!(truncate(_rotated_vec.z), 1.0);
    }

    #[test]
    fn rototranslation_apply_inverse_mut() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = 180.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(_origin_vector, _start_vector, _end_vector, _theta_rad);

        let mut _candidate_vector: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        let _rotated_vec = _roto.apply_inverse_mut(&mut _candidate_vector);

        assert_eq!(truncate(_rotated_vec.x), 3.0);
        assert_eq!(truncate(_rotated_vec.y), 2.0);
        assert_eq!(truncate(_rotated_vec.z), 1.0);
    }


    #[test]
    fn rotate_cube_around_axis_111() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(1.0, 1.0, 1.0);
        let angle = 2.0 * std::f64::consts::PI / 3.0;
        let rot = Rototranslation::around_axis(corner_1, corner_1, corner_2, angle);

        // println!("{:?}", rot.rotation_matrix());
        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        rot.apply_mut(&mut another_vec);
        println!("{}", another_vec);
        assert_eq!(truncate(another_vec.x), 0.0);
        assert_eq!(truncate(another_vec.y), 1.0);
        assert_eq!(truncate(another_vec.z), 1.0);
    }

    #[test]
    fn rotate_cube_around_axis_100() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(0.0, 0.0, 1.0);
        let angle = std::f64::consts::PI / 2.0;
        let rot = Rototranslation::around_axis(corner_1, corner_1, corner_2, angle);

        let mut another_vec = Vec3::new(1.0, 0.0, 0.0);
        rot.apply_mut(&mut another_vec);
        println!("{}", another_vec);
        assert_eq!(truncate(another_vec.x), 0.0);
        assert_eq!(truncate(another_vec.y), 1.0);
        assert_eq!(truncate(another_vec.z), 0.0);
    }
}
