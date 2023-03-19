extern crate bioshell_numerical;

#[cfg(test)]
mod matrix_test
{
    use bioshell_numerical::Vec3;
    use bioshell_numerical::matrix::Matrix3x3;

    #[test]
    fn matrix_from_column_values()
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
}