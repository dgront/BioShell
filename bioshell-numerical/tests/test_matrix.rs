//extern crate bioshell_numerical;

#[cfg(test)]
mod matrix_test {
    use bioshell_numerical::matrix::Matrix3x3;
    use bioshell_numerical::vec3::Vec3;

    #[test]
    fn matrix_unit_matrix() {
        let _vx = Vec3::new(1.0, 0.0, 0.0);
        let _vy = Vec3::new(0.0, 1.0, 0.0);
        let _vz = Vec3::new(1.0, 0.0, 1.0);
        let unit_mtx = Matrix3x3::from_column_vectors(&_vx, &_vy, &_vz);
        assert_eq!(unit_mtx[0], 1.0); assert_eq!(unit_mtx[4], 1.0); assert_eq!(unit_mtx[8], 1.0);
    }

    #[test]
    fn matrix_equality()
    {
        let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0]);
        assert_ne!(lhs, rhs);
        assert_eq!(lhs, lhs);
        assert_eq!(rhs, rhs);
    }

    #[test]
    fn matrix_add() {
        let mut lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        lhs.add(&rhs);
        let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_sub() {
        let mut lhs: Matrix3x3 =
            Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        lhs.sub(&rhs);
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_mul_scalar() {
        let mut lhs: Matrix3x3 =
            Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs = 2.0;
        lhs.mul_scalar(rhs);
        let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_div_scalar() {
        let mut lhs: Matrix3x3 =
            Matrix3x3::from_array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]);
        let rhs = 2.0;
        lhs.div_scalar(rhs);
        let expected = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_mul_vec_mut() {
        let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let mut rhs = Vec3::new(2.0, 2.0, 2.0);
        lhs.mul_vec_mut(&mut rhs);
        assert_eq!(12.0, rhs[0]);
        assert_eq!(30.0, rhs[1]);
        assert_eq!(48.0, rhs[2]);
    }

    #[test]
    fn matrix_mul_mat_mut() {
        let mut lhs: Matrix3x3 =
            Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        lhs.mul_mat_mut(&rhs);
        let expected = [30.0, 36.0, 42.0, 66.0, 81.0, 96.0, 102.0, 126.0, 150.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_det() {
        let mat: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let _det: f64 = mat.det();
        assert_eq!(0.0, _det);
    }

    #[test]
    fn matrix_identity()
    {
        let lhs = Matrix3x3::identity();
        let rhs = Matrix3x3::from_values(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn matrix_inverse() {
        let mut lhs: Matrix3x3 =
            Matrix3x3::from_array([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
        lhs.inverse();
        let expected = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(lhs[i], *val);
        }
    }

    #[test]
    fn matrix_from_values() {
        let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let mat: Matrix3x3 =
            Matrix3x3::from_values(cx.x, cx.y, cx.z, cy.x, cy.y, cy.z, cz.x, cz.y, cz.z);
        assert_eq!(cx.x, mat[0]); assert_eq!(cx.y, mat[1]); assert_eq!(cx.z, mat[2]);
        assert_eq!(cy.x, mat[3]); assert_eq!(cy.y, mat[4]); assert_eq!(cy.z, mat[5]);
        assert_eq!(cz.x, mat[6]); assert_eq!(cz.y, mat[7]); assert_eq!(cz.z, mat[8]);
    }

    #[test]
    fn matrix_from_column_vectors() {
        let _cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let _cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let _cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let _mat: Matrix3x3 = Matrix3x3::from_column_vectors(&_cx, &_cy, &_cx);
        assert_eq!(_cx.x, _mat[0]); assert_eq!(_cx.y, _mat[1]); assert_eq!(_cx.z, _mat[2]);
        assert_eq!(_cy.x, _mat[3]); assert_eq!(_cy.y, _mat[4]); assert_eq!(_cy.z, _mat[5]);
        assert_eq!(_cz.x, _mat[6]); assert_eq!(_cz.y, _mat[7]); assert_eq!(_cz.z, _mat[8]);
    }

    #[test]
    fn matrix_from_row_vectors() {
        let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let mat: Matrix3x3 = Matrix3x3::from_row_vectors(&cx, &cy, &cx);
        assert_eq!(cx.x, mat[0]); assert_eq!(cy.x, mat[1]); assert_eq!(cz.x, mat[2]);
        assert_eq!(cx.y, mat[3]); assert_eq!(cy.y, mat[4]); assert_eq!(cz.y, mat[5]);
        assert_eq!(cx.z, mat[6]); assert_eq!(cy.z, mat[7]); assert_eq!(cz.z, mat[8]);
    }

    #[test]
    fn matrix_add_s() {
        let _lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let _rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let _out = Matrix3x3::add_s(_lhs, _rhs);
        let _expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
        for (i, val) in _expected.iter().enumerate() {
            assert_eq!(_out[i], *val);
        }
    }

    #[test]
    fn matrix_sub_s() {
        let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let out = Matrix3x3::sub_s(lhs, rhs);
        let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(out[i], *val);
        }
    }

    #[test]
    fn matrix_mul_mat_s() {
        let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::from_array([9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
        let out = Matrix3x3::mul_mat_s(&lhs, &rhs);
        let expected = [30.0, 24.0, 18.0, 84.0, 69.0, 54.0, 138.0, 114.0, 90.0];
        for (i, val) in expected.iter().enumerate() {
            assert_eq!(out[i], *val);
        }
    }
}
