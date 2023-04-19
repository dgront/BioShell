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
        let _theta_rad = 360.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(_origin_vector, _start_vector, _end_vector, _theta_rad);

        let _candidate_vector: Vec3 = Vec3::new(3.0, 3.0, 3.0);

        let _rotated_vec = _roto.apply_mut(&_candidate_vector);

        assert_eq!(truncate(_rotated_vec.x), 3.0);
        assert_eq!(truncate(_rotated_vec.y), 3.0);
        assert_eq!(truncate(_rotated_vec.z), 3.0);
    }

    #[test]
    fn rototranslation_apply_inverse_mut() {
        let _origin_vector: Vec3 = Vec3::new(1.0, 1.0, 1.0);
        let _start_vector: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        let _end_vector: Vec3 = Vec3::new(4.0, 4.0, 4.0);
        let _theta_rad = 360.0 * std::f64::consts::PI /180.0;

        let _roto = Rototranslation::around_axis(_origin_vector, _start_vector, _end_vector, _theta_rad);

        let mut _candidate_vector: Vec3 = Vec3::new(3.0, 3.0, 3.0);

        let _rotated_vec = _roto.apply_inverse_mut(&mut _candidate_vector);

        assert_eq!(_rotated_vec.x, 3.0);
        assert_eq!(_rotated_vec.y, 3.0);
        assert_eq!(_rotated_vec.z, 3.0);
    }
}
