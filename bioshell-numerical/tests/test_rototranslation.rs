#[cfg(test)]
mod rototranslation_test {
    use bioshell_numerical::matrix::Matrix3x3;
    use bioshell_numerical::rototranslation::*;
    use bioshell_numerical::vec3::Vec3;

    #[test]
    fn rototranslation_struct() {
        let vx = Vec3::new(1.0, 0.0, 0.0);
        let vy = Vec3::new(0.0, 1.0, 0.0);
        let vz = Vec3::new(1.0, 0.0, 1.0);
        let unit_vec = Vec3::new(1.0, 1.0, 1.0);
        let another_vec = Vec3::new(10.0, 0.0, 10.0);
        let unit_mtx = Matrix3x3::from_values(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

        let rototran = Rototranslation::new(unit_mtx, unit_vec);

        rototran.apply(&another_vec);

        assert_eq!(10.0, another_vec[0]);
        assert_eq!(0.0, another_vec[1]);
        assert_eq!(10.0, another_vec[2]);
    }

    #[test]
    fn rotate_cube_around_axis_111() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(1.0, 1.0, 1.0);
        let angle = 2.0 * std::f64::consts::PI / 3.0;
        let rot = Rototranslation::around_axis(&corner_1, &corner_1, &corner_2, angle);

        println!("{:?}", rot.rotation_matrix());
        let mut another_vec = Vec3::new(1.0, 1.0, 0.0);
        for i in 0..4 {
            rot.apply_mut(&mut another_vec);
            println!("{}", another_vec);
        }


        // assert_eq!(10.0, another_vec.x);
        // assert_eq!(1.0, another_vec.y);
        // assert_eq!(10.0, another_vec.z);
    }

    #[test]
    fn rotate_cube_around_axis_001() {
        let corner_1 = Vec3::new(0.0, 0.0, 0.0);
        let corner_2 = Vec3::new(0.0, 0.0, 1.0);
        let angle = std::f64::consts::PI / 2.0;
        let rot = Rototranslation::around_axis(&corner_1, &corner_1, &corner_2, angle);

        let mut another_vec = Vec3::new(1.0, 1.0, 0.0);
        let  expected = vec![Vec3::new(-1.0, 1.0, 0.0), Vec3::new(-1.0, -1.0, 0.0),
                                 Vec3::new(1.0, -1.0, 0.0), Vec3::new(1.0, 1.0, 0.0)];

        for v in &expected {
            rot.apply_mut(&mut another_vec);
            assert!(v.distance_to(&another_vec).abs() < 0.0001);
        }

        let  expected_inv = vec![Vec3::new(1.0, -1.0, 0.0), Vec3::new(-1.0, -1.0, 0.0),
                             Vec3::new(-1.0, 1.0, 0.0), Vec3::new(1.0, 1.0, 0.0)];

        for v in &expected_inv {
            rot.apply_inverse_mut(&mut another_vec);
            assert!(v.distance_to(&another_vec).abs() < 0.0001);
        }
    }

    #[test]
    fn rototranslation_apply_mut() {
        let rotation_mat =
            Matrix3x3::from_values(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let translation_vec = Vec3::new(1.0, 1.0, 1.0);
        let mut another_vec = Vec3::new(10.0, 1.0, 10.0);

        let rot = Rototranslation::new(rotation_mat, translation_vec);

        rot.apply_mut(&mut another_vec);

        assert_eq!(10.0, another_vec.x);
        assert_eq!(1.0, another_vec.y);
        assert_eq!(10.0, another_vec.z);
    }

    #[test]
    fn rototranslation_apply() {
        let unit_matrix = Matrix3x3::from_values(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let vector = Vec3::new(10.0, 0.0, 10.0);

        let rot = Rototranslation::new(unit_matrix, vector);
        let rotated_vector = rot.apply(&vector);

        assert_eq!(rotated_vector.x, 10.0);
        assert_eq!(rotated_vector.y, 0.0);
        assert_eq!(rotated_vector.z, 10.0);
    }

    #[test]
    fn rototranslation_apply_inverse_mut() {
        let unit_matrix = Matrix3x3::from_values(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let trans_vec = Vec3::new(1.0, 1.0, 1.0);
        let rt = Rototranslation::new(unit_matrix, trans_vec);
        let mut another_vec = Vec3::new(1.0, 2.0, 3.0);

        rt.apply_inverse_mut(&mut another_vec);

        assert_eq!(Vec3::new(1.0, 2.0, 3.0), another_vec);
    }

    #[test]
    fn rototranslation_apply_inverse() {
        let unit_matrix = Matrix3x3::from_values(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        let trans_vec = Vec3::new(1.0, 1.0, 1.0);
        let rt = Rototranslation::new(unit_matrix, trans_vec);
        let another_vec = Vec3::new(1.0, 2.0, 3.0);

        let vec_out = rt.apply_inverse(&another_vec);

        assert_eq!(Vec3::new(1.0, 2.0, 3.0), vec_out);
    }
}
