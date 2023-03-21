//extern crate bioshell_numerical;

#[cfg(test)]
mod matrix_test
{
    use bioshell_numerical::matrix::Matrix3x3;
    use bioshell_numerical::vec3::Vec3;

    //region fn matrix_unit_matrix()
    #[test]
    fn matrix_unit_matrix()
    {
        let vx = Vec3::new(1.0, 0.0, 0.0);
        let vy = Vec3::new(0.0, 1.0, 0.0);
        let vz = Vec3::new(1.0, 0.0, 1.0);
        let unit_mtx = Matrix3x3::from_column_vectors(&vx, &vy, &vz);

        assert_eq!(unit_mtx[0], 1.0);
        assert_eq!(unit_mtx[4], 1.0);
        assert_eq!(unit_mtx[8], 1.0);
    }
    //endregion

    //region fn matrix_add()
    #[test]
    fn matrix_add()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        lhs.add(&rhs);

        assert_eq!(2.0, lhs[0]);
        assert_eq!(4.0, lhs[1]);
        assert_eq!(6.0, lhs[2]);
        assert_eq!(8.0, lhs[3]);
        assert_eq!(10.0, lhs[4]);
        assert_eq!(12.0, lhs[5]);
        assert_eq!(14.0, lhs[6]);
        assert_eq!(16.0, lhs[7]);
        assert_eq!(18.0, lhs[8]);
    }
    //endregion

    //region fn matrix_sub()
    #[test]
    fn matrix_sub()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);

        lhs.sub(&rhs);

        assert_eq!(0.0, lhs[0]);
        assert_eq!(0.0, lhs[1]);
        assert_eq!(0.0, lhs[2]);
        assert_eq!(0.0, lhs[3]);
        assert_eq!(0.0, lhs[4]);
        assert_eq!(0.0, lhs[5]);
        assert_eq!(0.0, lhs[6]);
        assert_eq!(0.0, lhs[7]);
        assert_eq!(0.0, lhs[8]);
    }
    //endregion

    //region fn matrix_mul_scalar()
    #[test]
    fn matrix_mul_scalar()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs = 2.0;

        lhs.mul_scalar(rhs);

        assert_eq!(2.0, lhs[0]);
        assert_eq!(4.0, lhs[1]);
        assert_eq!(6.0, lhs[2]);
        assert_eq!(8.0, lhs[3]);
        assert_eq!(10.0, lhs[4]);
        assert_eq!(12.0, lhs[5]);
        assert_eq!(14.0, lhs[6]);
        assert_eq!(16.0, lhs[7]);
        assert_eq!(18.0, lhs[8]);
    }
    //endregion

    //region fn matrix_div_scalar()
    #[test]
    fn matrix_div_scalar()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]);
        let rhs = 2.0;

        lhs.div_scalar(rhs);

        assert_eq!(1.0, lhs[0]);
        assert_eq!(2.0, lhs[1]);
        assert_eq!(3.0, lhs[2]);
        assert_eq!(4.0, lhs[3]);
        assert_eq!(5.0, lhs[4]);
        assert_eq!(6.0, lhs[5]);
        assert_eq!(7.0, lhs[6]);
        assert_eq!(8.0, lhs[7]);
        assert_eq!(9.0, lhs[8]);
    }
    //endregion

    //region fn matrix_mul_vec_mut()
    #[test]
    fn matrix_mul_vec_mut()
    {
        let lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0]);
        let mut rhs = Vec3::new(2.0,2.0,2.0);
        lhs.mul_vec_mut(&mut rhs);
        assert_eq!(12.0, rhs[0]);
        assert_eq!(30.0, rhs[1]);
        assert_eq!(48.0, rhs[2]);
    }
    //endregion

    //region fn matrix_mul_mat_mut()
    #[test]
    fn matrix_mul_mat_mut()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0]);
        let rhs = Matrix3x3::new_arr([1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0]);
        lhs.mul_mat_mut(&rhs);

        assert_eq!(30.0, lhs[0]);
        assert_eq!(36.0, lhs[1]);
        assert_eq!(42.0, lhs[2]);
        assert_eq!(66.0, lhs[3]);
        assert_eq!(81.0, lhs[4]);
        assert_eq!(96.0, lhs[5]);
        assert_eq!(102.0, lhs[6]);
        assert_eq!(126.0, lhs[7]);
        assert_eq!(150.0, lhs[8]);
    }
    //endregion

    //region matrix_det()
    #[test]
    fn matrix_det()
    {
        let mat: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0]);

        let _det: f64 = mat.det();

        assert_eq!(0.0, _det);
    }
    //endregion

    //region matrix_inverse()
    #[test]
    fn matrix_inverse()
    {
        let mut lhs: Matrix3x3 = Matrix3x3::new_arr([0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 0.0]);

        lhs.inverse();

        assert_eq!(0.0, lhs[0]);
        assert_eq!(0.0, lhs[1]);
        assert_eq!(1.0, lhs[2]);
        assert_eq!(1.0, lhs[3]);
        assert_eq!(0.0, lhs[4]);
        assert_eq!(0.0, lhs[5]);
        assert_eq!(0.0, lhs[6]);
        assert_eq!(1.0, lhs[7]);
        assert_eq!(0.0, lhs[8]);
    }
    //endregion

    //region fn matrix_from_column_values()
    #[test]
    fn matrix_from_values()
    {
        let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        let mat: Matrix3x3 = Matrix3x3::from_values(cx.x, cx.y, cx.z,
                                                    cy.x, cy.y, cy.z,
                                                    cz.x, cz.y, cz.z);

        assert_eq!(cx.x, mat[0]);
        assert_eq!(cx.y, mat[1]);
        assert_eq!(cx.z, mat[2]);

        assert_eq!(cy.x, mat[3]);
        assert_eq!(cy.y, mat[4]);
        assert_eq!(cy.z, mat[5]);

        assert_eq!(cz.x, mat[6]);
        assert_eq!(cz.y, mat[7]);
        assert_eq!(cz.z, mat[8]);
    }
    //endregion

    //region fn matrix_from_column_vectors()
    #[test]
    fn matrix_from_column_vectors()
    {
        let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        let mat: Matrix3x3 = Matrix3x3::from_column_vectors(&cx, &cy, &cx);

        assert_eq!(cx.x, mat[0]);
        assert_eq!(cx.y, mat[1]);
        assert_eq!(cx.z, mat[2]);

        assert_eq!(cy.x, mat[3]);
        assert_eq!(cy.y, mat[4]);
        assert_eq!(cy.z, mat[5]);

        assert_eq!(cz.x, mat[6]);
        assert_eq!(cz.y, mat[7]);
        assert_eq!(cz.z, mat[8]);
    }
    //endregion

    //region fn matrix_from_row_vectors()
    #[test]
    fn matrix_from_row_vectors()
    {
        let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
        let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);

        let mat: Matrix3x3 = Matrix3x3::from_row_vectors(&cx, &cy, &cx);

        assert_eq!(cx.x, mat[0]);
        assert_eq!(cy.x, mat[1]);
        assert_eq!(cz.x, mat[2]);

        assert_eq!(cx.y, mat[3]);
        assert_eq!(cy.y, mat[4]);
        assert_eq!(cz.y, mat[5]);

        assert_eq!(cx.z, mat[6]);
        assert_eq!(cy.z, mat[7]);
        assert_eq!(cz.z, mat[8]);
    }
    //endregion

    //region fn matrix_add_s()
    #[test]
    fn matrix_add_s()
    {
        let lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let out= Matrix3x3::add_s(lhs, rhs);

        assert_eq!(2.0, out[0]);
        assert_eq!(4.0, out[1]);
        assert_eq!(6.0, out[2]);
        assert_eq!(8.0, out[3]);
        assert_eq!(10.0, out[4]);
        assert_eq!(12.0, out[5]);
        assert_eq!(14.0, out[6]);
        assert_eq!(16.0, out[7]);
        assert_eq!(18.0, out[8]);
    }
    //endregion

    //region fn matrix_sub_s()
    #[test]
    fn matrix_sub_s()
    {
        let lhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let rhs: Matrix3x3 = Matrix3x3::new_arr([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let out= Matrix3x3::sub_s(lhs, rhs);

        assert_eq!(0.0, out[0]);
        assert_eq!(0.0, out[1]);
        assert_eq!(0.0, out[2]);
        assert_eq!(0.0, out[3]);
        assert_eq!(0.0, out[4]);
        assert_eq!(0.0, out[5]);
        assert_eq!(0.0, out[6]);
        assert_eq!(0.0, out[7]);
        assert_eq!(0.0, out[8]);
    }
    //endregion
}