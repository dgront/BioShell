use std::fmt;
use std::ops::{Index, IndexMut, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::calc::Vec3;

/// Represents a 3x3 matrix, e.g. for linear 3D transformations.
///
/// Internally the matrix elements are stored as an `[f64; 9]` array in the order as follows:
/// ```text
///     | 0 1 2 |
/// m = | 3 4 5 |
///     | 6 7 8 |
/// ```
///
/// # Example
/// ```rust
/// use bioshell_pdb::calc::Matrix3x3;
/// use bioshell_pdb::calc::Vec3;
///
/// let vx = Vec3::new(0.0, 3.0, 6.0);
/// let vy = Vec3::new(1.0, 4.0, 7.0);
/// let vz = Vec3::new(2.0, 5.0, 8.0);
/// let unit_mtx = Matrix3x3::from_column_vectors(&vx, &vy, &vz);
/// assert_eq!(unit_mtx[0], 0.0); assert_eq!(unit_mtx[3], 3.0); assert_eq!(unit_mtx[7], 7.0);
/// ```

#[derive(Clone, Copy, Default)]
pub struct Matrix3x3 {
    array: [f64; 9],
}

impl Index<usize> for Matrix3x3 {
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64 {
        &self.array[i]
    }
}

impl IndexMut<usize> for Matrix3x3 {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut f64 {
        &mut self.array[i]
    }
}

impl fmt::Debug for Matrix3x3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f, "[ [{:.2}, {:.2}, {:.2}], [{:.2}, {:.2}, {:.2}], [{:.2}, {:.2}, {:.2}] ]",
            self.array[0], self.array[1], self.array[2], self.array[3], self.array[4],
            self.array[5], self.array[6], self.array[7], self.array[8]
        )
    }
}

impl AddAssign<&Matrix3x3> for Matrix3x3 {
    /// Provides `+=` operator that adds another matrix to this matrix.
    ///
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// lhs += &rhs;
    /// # let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
    /// # for (i, val) in expected.iter().enumerate() {
    /// #     assert_eq!(lhs[i], *val);
    /// # }
    /// ```
    fn add_assign(&mut self, rhs: &Matrix3x3) {
        for i in 0..9 { self[i] += rhs[i]; }
    }
}

impl SubAssign<&Matrix3x3> for Matrix3x3 {
    /// Provides `-=` operator that subtracts another matrix from this matrix.
    ///
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// lhs -= &rhs;
    /// # for i in 0..9 {
    /// #     assert_eq!(lhs[i], 0.0);
    /// # }
    /// ```
    fn sub_assign(&mut self, rhs: &Matrix3x3) {
        for i in 0..9 { self[i] -= rhs[i]; }
    }
}

impl DivAssign<f64> for Matrix3x3 {
    /// Provides `/=` operator that divides this matrix by a scalar value
    ///
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::identity();
    /// lhs /= 10.0;
    /// assert_eq!(lhs[4], 0.1);
    /// ```
    fn div_assign(&mut self, scalar: f64) {
        for i in 0..9 { self[i] /= scalar; }
    }
}

impl MulAssign<f64> for Matrix3x3 {
    /// Provides `*=` operator that multiplies this matrix by a scalar value
    ///
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::identity();
    /// lhs *= 10.0;
    /// assert_eq!(lhs[4], 10.0);
    /// ```
    fn mul_assign(&mut self, scalar: f64) {
        for i in 0..9 { self[i] *= scalar; }
    }
}

impl PartialEq for Matrix3x3 {
    fn eq(&self, other: &Self) -> bool {
        let mut returns: bool = true;

        for i in 1..9 {
            returns = returns && (self[i] == other[i]);
        }

        return returns;
    }
}

impl Matrix3x3 {
    /// Creates and returns a new Matrix3x3 instance with all elements initialized to 0.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    /// let mat = Matrix3x3::new();
    /// assert_eq!(Matrix3x3::from_array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,]), mat)
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Constructs a new Matrix3x3 object from an array of 9 f64 elements.
    ///
    /// # Arguments
    ///
    /// * m - An array of 9 `f64` elements, which are used to fill the matrix in column-major order.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let m = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    /// let mtx = Matrix3x3::from_array(m);
    ///
    /// assert_eq!(mtx[0], 1.0);
    /// assert_eq!(mtx[1], 2.0);
    /// assert_eq!(mtx[2], 3.0);
    /// ```
    pub fn from_array(m: [f64; 9]) -> Self {
        Matrix3x3 { array: m }
    }

    /// Creates a Matrix3x3 from three column vectors, representing the columns of the matrix.
    ///
    /// Given three column vectors ``a``, ``b`` and ``c``, this method returns a matrix created as:
    /// ```math
    /// M = \begin{vmatrix}
    /// a & b & c\\
    /// \end{vmatrix} = \begin{vmatrix}
    /// a_x & b_x & c_x\\
    /// a_y & b_y & c_y\\
    /// a_z & b_z & c_z\\
    /// \end{vmatrix}
    /// ```
    ///
    /// # Arguments
    /// * a - a Vec3 object representing the first column of the matrix.
    /// * b - a Vec3 object representing the second column of the matrix.
    /// * c - a Vec3 object representing the third column of the matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::{Vec3, Matrix3x3};
    ///
    /// let a: Vec3 = Vec3::new(1.0, 4.0, 7.0);
    /// let b: Vec3 = Vec3::new(2.0, 5.0, 8.0);
    /// let c: Vec3 = Vec3::new(3.0, 6.0, 9.0);
    /// let mat: Matrix3x3 = Matrix3x3::from_column_vectors(&a, &b, &c);
    /// assert_eq!(a.x, mat[0]); assert_eq!(a.y, mat[3]); assert_eq!(a.z, mat[6]);
    /// assert_eq!(b.x, mat[1]); assert_eq!(b.y, mat[4]); assert_eq!(b.z, mat[7]);
    /// assert_eq!(c.x, mat[2]); assert_eq!(c.y, mat[5]); assert_eq!(c.z, mat[8]);
    ///```
    pub fn from_column_vectors(a: &Vec3, b: &Vec3, c: &Vec3) -> Self {
        Self::from_array([a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z])
    }

    /// Creates a Matrix3x3 from three column vectors, representing the columns of the matrix.
    ///
    /// Given three row vectors ``a``, ``b`` and ``c``, this method returns a matrix created as:
    /// ```math
    /// M = \begin{vmatrix}
    /// a \\
    /// b \\
    /// c \\
    /// \end{vmatrix} = \begin{vmatrix}
    /// a_x & a_y & a_z\\
    /// b_x & b_y & b_z\\
    /// c_x & c_y & c_z\\
    /// \end{vmatrix}
    /// ```
    ///
    /// # Arguments
    /// * a - a Vec3 object representing the first row of the matrix.
    /// * b - a Vec3 object representing the second row of the matrix.
    /// * c - a Vec3 object representing the third column of the matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::{Vec3, Matrix3x3};
    ///
    /// let a: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let b: Vec3 = Vec3::new(4.0, 5.0, 6.0);
    /// let c: Vec3 = Vec3::new(7.0, 8.0, 9.0);
    /// let mat: Matrix3x3 = Matrix3x3::from_row_vectors(&a, &b, &c);
    /// assert_eq!(a.x, mat[0]); assert_eq!(a.y, mat[1]); assert_eq!(a.z, mat[2]);
    /// assert_eq!(b.x, mat[3]); assert_eq!(b.y, mat[4]); assert_eq!(b.z, mat[5]);
    /// assert_eq!(c.x, mat[6]); assert_eq!(c.y, mat[7]); assert_eq!(c.z, mat[8]);
    ///```
    pub fn from_row_vectors(a: &Vec3, b: &Vec3, c: &Vec3) -> Self {
        Self::from_array([a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z])
    }

    /// Provides access to an elements of this matrix for a given ``row`` and ``column`` indexes.
    ///
    /// Effectively this method provides a 2D indexing capability for this 1D data storage.
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::assert_delta;
    /// use bioshell_pdb::calc::Matrix3x3;
    /// let m = Matrix3x3::from_array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]);
    /// assert_delta!(m.elem(0, 1), 1.0, 0.0000001);
    /// assert_delta!(m.elem(0, 1), 1.0, 0.0000001);
    /// assert_delta!(m.elem(1, 2), 5.0, 0.0000001);
    /// ```
    pub fn elem(&self, row: usize, col:usize) -> f64 {
        self[row * 3 + col]
    }

    /// Multiplies the caller matrix object with a Vec3 in place, i.e., modifies the matrix object itself.
    ///
    /// `rhs` - is a mutable reference to a Vec3 object that represents the vector to be multiplied with the matrix.
    /// The multiplication is performed as M x v, where M is the matrix object and v is the Vec3 object.
    /// The result of this operation is stored back in rhs.
    ///
    /// # Arguments
    ///
    /// * `rhs` - is a mutable reference to a Vec3 object that represents the vector to be multiplied with the matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::assert_vec3_eq;
    /// use bioshell_pdb::calc::{Vec3, Matrix3x3};
    ///
    /// let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let mut rhs = Vec3::new(2.0, 2.0, 2.0);
    /// lhs.mul_vec_mut(&mut rhs);
    /// assert_vec3_eq!(Vec3::new(12.0, 30.0, 48.0), rhs, 0.0001, "Incorrect vector after multiplication");
    ///```
    pub fn mul_vec_mut(&self, rhs: &mut Vec3) {
        let lhs = self;
        let x = lhs[0] * rhs.x + lhs[1] * rhs.y + lhs[2] * rhs.z;
        let y = lhs[3] * rhs.x + lhs[4] * rhs.y + lhs[5] * rhs.z;
        let z = lhs[6] * rhs.x + lhs[7] * rhs.y + lhs[8] * rhs.z;
        rhs.x = x;
        rhs.y = y;
        rhs.z = z;
    }

    /// Multiplies the caller matrix object with another matrix in-place.
    ///
    /// The multiplication operation takes the current matrix object (left-hand-side) and the given matrix (right-hand-side),
    /// and computes their matrix product, i.e. the resulting matrix has the same number of rows as the left-hand-side matrix,
    /// and the same number of columns as the right-hand-side matrix.
    ///
    /// Each element of the resulting matrix is obtained by taking the dot product of the corresponding row of the left-hand-side
    /// matrix and the corresponding column of the right-hand-side matrix.
    ///
    /// # Arguments
    ///
    /// * rhs - a reference to the matrix object to be multiplied with the caller matrix object.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::{Matrix3x3};
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// lhs.mul_mat_mut(&rhs);
    /// let expected = [30.0, 36.0, 42.0, 66.0, 81.0, 96.0, 102.0, 126.0, 150.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    /// ```
    pub fn mul_mat_mut(&mut self, rhs: &Matrix3x3) {
        let lhs = self;
        let temp = [
            lhs[0] * rhs[0] + lhs[1] * rhs[3] + lhs[2] * rhs[6], // row 1, column 1
            lhs[0] * rhs[1] + lhs[1] * rhs[4] + lhs[2] * rhs[7], // row 1, column 2
            lhs[0] * rhs[2] + lhs[1] * rhs[5] + lhs[2] * rhs[8], // row 1, column 3
            lhs[3] * rhs[0] + lhs[4] * rhs[3] + lhs[5] * rhs[6], // row 2, column 1
            lhs[3] * rhs[1] + lhs[4] * rhs[4] + lhs[5] * rhs[7], // row 2, column 2
            lhs[3] * rhs[2] + lhs[4] * rhs[5] + lhs[5] * rhs[8], // row 2, column 3
            lhs[6] * rhs[0] + lhs[7] * rhs[3] + lhs[8] * rhs[6], // row 3, column 1
            lhs[6] * rhs[1] + lhs[7] * rhs[4] + lhs[8] * rhs[7], // row 3, column 2
            lhs[6] * rhs[2] + lhs[7] * rhs[5] + lhs[8] * rhs[8], // row 3, column 3
        ];

        for i in 0..9 {
            lhs[i] = temp[i];
        }
    }

    /// Computes the determinant of the 3x3 matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::{Matrix3x3};
    ///
    /// let mat: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let _det: f64 = mat.det();
    /// assert_eq!(0.0, _det);
    ///```
    pub fn det(&self) -> f64 {
        let lhs = &self.array;
        let a = lhs[0] * (lhs[4] * lhs[8] - lhs[5] * lhs[7]);
        let b = lhs[1] * (lhs[3] * lhs[8] - lhs[5] * lhs[6]);
        let c = lhs[2] * (lhs[3] * lhs[7] - lhs[4] * lhs[6]);
        return a - b + c;
    }

    /// Calculates the inverse of the matrix in-place.
    ///
    /// # Panics
    ///
    /// This function may panic if the matrix is not invertible, i.e. has a determinant of zero.
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::from_array([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
    /// lhs.inverse();
    /// let expected = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    ///```
    /// todo: Throw exception when det==0
    pub fn inverse(&mut self) {
        let _det = self.det();
        let lhs = self;
        if _det != 0.0 {
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

            for i in 0..9 {
                lhs[i] = temp[i] * inv_det;
            }
        }
    }

    /// Transposes this matrix in-place.
    pub fn transpose(&mut self) {
        let lhs = self;
        let temp = [lhs[0], lhs[3], lhs[6], lhs[1], lhs[4], lhs[7], lhs[2], lhs[5], lhs[8], ];
        for i in 0..9 { lhs[i] = temp[i]; }
    }

    /// Adds two Matrix3x3 matrices and returns the result.
    ///
    /// # Arguments
    ///
    /// * lhs - A Matrix3x3 instance representing the left-hand side of the operation.
    /// * rhs - A Matrix3x3 instance representing the right-hand side of the operation.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let out = Matrix3x3::add_s(lhs, rhs);
    /// let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(out[i], *val);
    /// }
    /// ```
    pub fn add_s(lhs: Matrix3x3, rhs: Matrix3x3) -> Matrix3x3
    {
        return Matrix3x3::from_array([
            lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],
            lhs[3] + rhs[3], lhs[4] + rhs[4], lhs[5] + rhs[5],
            lhs[6] + rhs[6], lhs[7] + rhs[7], lhs[8] + rhs[8],
        ]);
    }

    /// Performs matrix subtraction of two Matrix3x3 objects.
    ///
    /// Subtracts the matrix rhs from the matrix lhs and returns the resulting Matrix3x3.
    ///
    /// # Arguments
    ///
    /// * lhs - The matrix object from which rhs is to be subtracted.
    /// * rhs - The matrix object which is to be subtracted from lhs.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let out = Matrix3x3::sub_s(lhs, rhs);
    /// let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(out[i], *val);
    /// }
    ///```
    pub fn sub_s(lhs: Matrix3x3, rhs: Matrix3x3) -> Matrix3x3 {
        return Matrix3x3::from_array([
            lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2],
            lhs[3] - rhs[3], lhs[4] - rhs[4], lhs[5] - rhs[5],
            lhs[6] - rhs[6], lhs[7] - rhs[7], lhs[8] - rhs[8],
        ]);
    }

    /// Returns a 3x3 identity matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let lhs = Matrix3x3::identity();
    /// let rhs = Matrix3x3::from_array([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    /// assert_eq!(lhs, rhs);
    /// ```
    pub fn identity() -> Self {
        return Matrix3x3::from_array(  [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
    }

    /// Multiplies a 3x3 matrix with a 3D vector and returns the resulting 3D vector.
    ///
    /// Returns a newly created vector that is the result of a matrix-by-vector multiplication
    /// ```text
    /// |   |   |     |   |   |
    /// | v | = |  M  | * | v |
    /// |   |   |     |   |   |
    /// ```
    /// # Arguments
    ///
    /// * mat - A reference to a Matrix3x3 object.
    /// * vec - A reference to a Vec3 object.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::{Vec3, Matrix3x3};
    ///
    /// let mat = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let vec = Vec3::new(1.0, 2.0, 3.0);
    /// let result = Matrix3x3::mul_vec_s(&mat, &vec);
    /// assert_eq!(result, Vec3::new(14.0, 32.0, 50.0));
    ///```
    pub fn mul_vec_s(mat: &Matrix3x3, vec: &Vec3) -> Vec3 {
        let mut vec_copy = vec.clone();
        mat.mul_vec_mut(&mut vec_copy);
        return vec_copy;
    }

    /// Multiplies two Matrix3x3 objects and returns the result as a new Matrix3x3.
    ///
    /// This function does not modify the input matrices, but rather creates and returns a new Matrix3x3 object
    /// which is the result of multiplying the lhs matrix with the rhs matrix.
    ///
    /// # Arguments
    ///
    /// * lhs - A reference to the left-hand-side Matrix3x3 object.
    /// * rhs - A reference to the right-hand-side Matrix3x3 object.
    ///
    /// # Returns
    ///
    /// A new Matrix3x3 object which is the result of multiplying the lhs matrix with the rhs matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_pdb::calc::Matrix3x3;
    ///
    /// let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]);
    /// let out = Matrix3x3::mul_mat_s(&lhs, &rhs);
    /// let expected = [30.0, 24.0, 18.0, 84.0, 69.0, 54.0, 138.0, 114.0, 90.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(out[i], *val);
    /// }
    ///```
    pub fn mul_mat_s(lhs: &Matrix3x3, rhs: &Matrix3x3) -> Matrix3x3 {
        let mut lhs_copy = lhs.clone();
        lhs_copy.mul_mat_mut(&rhs);
        return lhs_copy;
    }

    pub fn outer_product_mat(a: &Matrix3x3, b: &Matrix3x3) -> Matrix3x3 {
        let mut result = Matrix3x3::default();
        for i in 0..3 {
            for j in 0..3 {
                result[i * 3 + j] = a[i] * b[j];
            }
        }
        return result;
    }
}
