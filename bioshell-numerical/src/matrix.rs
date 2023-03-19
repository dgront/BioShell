use std::ops::{Index, IndexMut};
use crate::vec3::Vec3;

/// Represents a 3x3 matrix, e.g. for linear transformations.
///
/// Internally the matrix elements are stored as an `[f64; 9]` array.
///
/// # Example
/// ```rust
/// use bioshell_numerical::matrix::Matrix3x3;
/// use bioshell_numerical::Vec3;
/// let vx = Vec3::new(1.0, 0.0, 0.0);
/// let vy = Vec3::new(0.0, 1.0, 0.0);
/// let vz = Vec3::new(0.0, 0.0, 1.0);
/// let unit_mtx = Matrix3x3::from_column_vectors(&vx, &vy, &vz);
///
/// assert_eq!(unit_mtx[0], 1.0);
/// assert_eq!(unit_mtx[4], 1.0);
/// assert_eq!(unit_mtx[8], 1.0);
/// ```
//region pub struct Matrix3x3
#[derive(Clone)]
pub struct Matrix3x3
{
    _array: [f64; 9],
}
//endregion

//region Indexing facility implementation.
impl Index<usize> for Matrix3x3
{
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64
    {
        &self._array[i]
    }
}

impl IndexMut<usize> for Matrix3x3
{
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut f64
    {
        &mut self._array[i]
    }
}
//endregion

impl Matrix3x3
{
    // public member functions
    //region pub fn show(&self)
    //display the matrix
    pub fn show(&self)
    {
        println!("[[{} {} {}]", self[0], self[1], self[2]);
        println!(" [{} {} {}]", self[3], self[4], self[5]);
        println!(" [{} {} {}]]", self[6], self[7], self[8]);
    }
    //endregion

    //region pub fn add(&self, rhs: Matrix3x3)
    //Add this matrix to another matrix.
    pub fn add(&mut self, rhs: &Matrix3x3)
    {
        let lhs = self;
        lhs[0] = lhs[0] + rhs[0];
        lhs[1] = lhs[1] + rhs[1];
        lhs[2] = lhs[2] + rhs[2];
        lhs[3] = lhs[3] + rhs[3];
        lhs[4] = lhs[4] + rhs[4];
        lhs[5] = lhs[5] + rhs[5];
        lhs[6] = lhs[6] + rhs[6];
        lhs[7] = lhs[7] + rhs[7];
        lhs[8] = lhs[8] + rhs[8];
    }
    //endregion

    //region pub fn sub(&self, rhs: Matrix3x3)
    //subtract the RHS matrix from the caller matrix object.
    pub fn sub(&mut self, rhs: Matrix3x3)
    {
        let lhs = self;
        lhs[0] = lhs[0] - rhs[0];
        lhs[1] = lhs[1] - rhs[1];
        lhs[2] = lhs[2] - rhs[2];
        lhs[3] = lhs[3] - rhs[3];
        lhs[4] = lhs[4] - rhs[4];
        lhs[5] = lhs[5] - rhs[5];
        lhs[6] = lhs[6] - rhs[6];
        lhs[7] = lhs[7] - rhs[7];
        lhs[8] = lhs[8] - rhs[8];
    }
    //endregion

    //region pub fn mul_scalar(&self, rhs:f64)
    //multiply the caller matrix object with a scalar value.
    pub fn mul_scalar(&mut self, rhs:f64)
    {
        let lhs = self;
        lhs[0] = lhs[0] * rhs;
        lhs[1] = lhs[1] * rhs;
        lhs[2] = lhs[2] * rhs;
        lhs[3] = lhs[3] * rhs;
        lhs[4] = lhs[4] * rhs;
        lhs[5] = lhs[5] * rhs;
        lhs[6] = lhs[6] * rhs;
        lhs[7] = lhs[7] * rhs;
        lhs[8] = lhs[8] * rhs;
    }
    //endregion

    //region pub fn div_scalar(&self, rhs:f64)
    //divide the caller matrix object with a scalar value.
    pub fn div_scalar(&mut self, rhs:f64)
    {
        let lhs = self;
        lhs[0] = lhs[0] / rhs;
        lhs[1] = lhs[1] / rhs;
        lhs[2] = lhs[2] / rhs;
        lhs[3] = lhs[3] / rhs;
        lhs[4] = lhs[4] / rhs;
        lhs[5] = lhs[5] / rhs;
        lhs[6] = lhs[6] / rhs;
        lhs[7] = lhs[7] / rhs;
        lhs[8] = lhs[8] / rhs;
    }
    //endregion

    //region pub fn mul_vec_mut(&self, v: &mut Vec3)
    //Multiply the caller matrix object with a vector object.
    pub fn mul_vec_mut(&self, v: &mut Vec3)
    {
        let lhs = self;
        let x = lhs[0] * v.x + lhs[1] * v.y + lhs[2] * v.z;
        let y = lhs[3] * v.x + lhs[4] * v.y + lhs[5] * v.z;
        let z = lhs[6] * v.x + lhs[7] * v.y + lhs[8] * v.z;
        v.x = x;
        v.y = y;
        v.z = z;
    }
    //endregion

    //region pub fn mul_mat_mut(&self, rhs: &mut Matrix3x3)
    //Multiply the caller matrix object with another matrix object.
    //The result is returned to the RHS.
    pub fn mul_mat_mut(&self, rhs: &mut Matrix3x3)
    {
        let lhs = self;
        let temp = [
            lhs[0] * rhs[0] + lhs[1] * rhs[3] + lhs[2] * rhs[6],  // row 1, column 1
            lhs[0] * rhs[1] + lhs[1] * rhs[4] + lhs[2] * rhs[7],  // row 1, column 2
            lhs[0] * rhs[2] + lhs[1] * rhs[5] + lhs[2] * rhs[8],  // row 1, column 3
            lhs[3] * rhs[0] + lhs[4] * rhs[3] + lhs[5] * rhs[6],  // row 2, column 1
            lhs[3] * rhs[1] + lhs[4] * rhs[4] + lhs[5] * rhs[7],  // row 2, column 2
            lhs[3] * rhs[2] + lhs[4] * rhs[5] + lhs[5] * rhs[8],  // row 2, column 3
            lhs[6] * rhs[0] + lhs[7] * rhs[3] + lhs[8] * rhs[6],  // row 3, column 1
            lhs[6] * rhs[1] + lhs[7] * rhs[4] + lhs[8] * rhs[7],  // row 3, column 2
            lhs[6] * rhs[2] + lhs[7] * rhs[5] + lhs[8] * rhs[8],  // row 3, column 3
        ];
        rhs[0] = temp[0];
        rhs[1] = temp[1];
        rhs[2] = temp[2];
        rhs[3] = temp[3];
        rhs[4] = temp[4];
        rhs[5] = temp[5];
        rhs[6] = temp[6];
        rhs[7] = temp[7];
        rhs[8] = temp[8];
    }
    //endregion

    //region pub fn det(&self) -> f64
    //Calculate the determinant value of a matrix.
    pub fn det(&self) -> f64
    {
        let lhs = self;
        return lhs[0] * (lhs[4] * lhs[8] - lhs[5] * lhs[7])
            - lhs[1] * (lhs[3] * lhs[8] - lhs[5] * lhs[6])
            + lhs[2] * (lhs[3] * lhs[7] - lhs[4] * lhs[6]);
    }
    //endregion

    //region pub fn inverse(&self)
    //Inverse the caller matrix object.
    pub fn inverse(&mut self)
    {
        let _det = self.det();
        let lhs = self;
        if _det != 0.0
        {
            let inv_det = 1.0 / _det;
            let temp = [
                lhs[4] * lhs[8] - lhs[5] * lhs[7],
                lhs[2] * lhs[7] - lhs[1] * lhs[8],
                lhs[1] * lhs[5] - lhs[2] * lhs[4],
                lhs[5] * lhs[6] - lhs[3] * lhs[8],
                lhs[0] * lhs[8] - lhs[2] * lhs[6],
                lhs[2] * lhs[3] - lhs[0] * lhs[5],
                lhs[3] * lhs[7] - lhs[4] * lhs[6],
                lhs[1] * lhs[6] - lhs[0] * lhs[7],
                lhs[0] * lhs[4] - lhs[1] * lhs[3],
            ];
            lhs[0] = temp[0] * inv_det;
            lhs[1] = temp[1] * inv_det;
            lhs[2] = temp[2] * inv_det;
            lhs[3] = temp[3] * inv_det;
            lhs[4] = temp[4] * inv_det;
            lhs[5] = temp[5] * inv_det;
            lhs[6] = temp[6] * inv_det;
            lhs[7] = temp[7] * inv_det;
            lhs[8] = temp[8] * inv_det;
        }
    }
    //endregion

    // static functions
    //region pub fn new(m: [f64; 9]) -> Self
    //Construct a matrix object from an array of scalar values.
    pub fn new(m: [f64; 9]) -> Self
    {
        Matrix3x3 { _array: m }
    }
    //endregion

    //region pub fn from_values()
    //Construct a matrix object from nine scalar values.
    pub fn from_values(cx_x:f64, cx_y:f64, cx_z:f64,
                       cy_x:f64, cy_y:f64, cy_z:f64,
                       cz_x:f64, cz_y:f64, cz_z:f64) ->Self
    {
        Self::new([
            cx_x, cx_y, cx_z,
            cy_x, cy_y, cy_z,
            cz_x, cz_y, cz_z,
        ])
    }
    //endregion

    //region pub fn from_column_vectors(cx: &Vec3, cy: &Vec3, cz: &Vec3) -> Self
    //Construct a matrix object from three column vector objects.
    pub fn from_column_vectors(cx: &Vec3, cy: &Vec3, cz: &Vec3) -> Self
    {
        let cx1 = cx.normalized();
        let cy1 = cy.normalized();
        let cz1 = cz.normalized();

        Self::new([
            cx1.x, cx1.y, cx1.z,
            cy1.x, cy1.y, cy1.z,
            cz1.x, cz1.y, cz1.z,
        ])
    }
    //endregion

    //region pub fn from_row_vectors(rx: &Vec3, ry: &Vec3, rz: &Vec3) -> Self
    //Construct a matrix object from three row vector objects.
    pub fn from_row_vectors(rx: &Vec3, ry: &Vec3, rz: &Vec3) -> Self
    {
        let rx1 = rx.normalized();
        let ry1 = ry.normalized();
        let rz1 = rz.normalized();

        Self::new([
            rx1.x, ry1.x, rz1.x,
            rx1.y, ry1.y, rz1.y,
            rx1.z, ry1.z, rz1.z,
        ])
    }
    //endregion

    //region pub fn add_s(lhs : Matrix3x3, rhs: Matrix3x3)-> Matrix3x3
    //Add two matrices and return the result as a new matrix.
    pub fn add_s(lhs : Matrix3x3, rhs: Matrix3x3) -> Matrix3x3
    {
        return Matrix3x3::new([lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3],
            lhs[4] + rhs[4], lhs[5] + rhs[5], lhs[6] + rhs[6], lhs[7] + rhs[7], lhs[8] + rhs[8]])
    }
    //endregion

    //region pub fn sub_s(lhs : Matrix3x3, rhs: Matrix3x3)
    //Subtract one matrix from another and return result as a new matrix.
    pub fn sub_s(lhs : Matrix3x3, rhs: Matrix3x3) -> Matrix3x3
    {
        return Matrix3x3::new([lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3],
            lhs[4] - rhs[4], lhs[5] - rhs[5], lhs[6] - rhs[6], lhs[7] - rhs[7], lhs[8] - rhs[8]])
    }
    //endregion

    //region pub fn inverse_s(mat: Matrix3x3)->Matrix3x3
    //Inverse a matrix and return the result as a new matrix.
    pub fn inverse_s(mat: Matrix3x3)->Matrix3x3
    {
        let mut mat_copy = mat.clone();
        mat_copy.inverse();
        return mat_copy;
    }
    //endregion

    //region pub fn mul_vec_s(mat: &Matrix3x3, vec: &Vec3)-> Vec3
    //Multiply a matrix with a vector and return the result as a new vector.
    pub fn mul_vec_s(mat: &Matrix3x3, vec: &Vec3)-> Vec3
    {
        let mut vec_copy = vec.clone();
        mat.mul_vec_mut(&mut vec_copy);
        return vec_copy;
    }
    //endregion

    //region pub fn mul_mat_s(lhs: &Matrix3x3, rhs:&Matrix3x3) -> Matrix3x3
    //Multiply two matrices and return the result as a new matrix.
    pub fn mul_mat_s(lhs: &Matrix3x3, rhs:&Matrix3x3) -> Matrix3x3
    {
        let mut rhs_copy = rhs.clone();
        lhs.mul_mat_mut(&mut rhs_copy);
        return rhs_copy;
    }
    //endregion
}