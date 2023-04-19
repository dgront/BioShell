use crate::vec3::Vec3;
use std::fmt;
use std::ops::{Index, IndexMut};
use std::ops::Mul;

/// Represents a 3x3 matrix, e.g. for linear transformations.
///
/// Internally the matrix elements are stored as an `[f64; 9]` array.
///
/// # Example
/// ```rust
/// use bioshell_numerical::matrix::Matrix3x3;
/// use bioshell_numerical::Vec3;
///
/// let vx = Vec3::new(1.0, 0.0, 0.0);
/// let vy = Vec3::new(0.0, 1.0, 0.0);
/// let vz = Vec3::new(1.0, 0.0, 1.0);
/// let unit_mtx = Matrix3x3::from_column_vectors(&vx, &vy, &vz);
/// assert_eq!(unit_mtx[0], 1.0); assert_eq!(unit_mtx[4], 1.0); assert_eq!(unit_mtx[8], 1.0);
/// ```

#[derive(Clone, Copy, Default)]
pub struct Matrix3x3 {
    _array: [f64; 9],
}

impl Index<usize> for Matrix3x3 {
    type Output = f64;
    fn index<'a>(&'a self, i: usize) -> &'a f64 {
        &self._array[i]
    }
}

impl IndexMut<usize> for Matrix3x3 {
    fn index_mut<'a>(&'a mut self, i: usize) -> &'a mut f64 {
        &mut self._array[i]
    }
}

impl fmt::Debug for Matrix3x3 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "[ [{:.2}, {:.2}, {:.2}], [{:.2}, {:.2}, {:.2}], [{:.2}, {:.2}, {:.2}] ]",
            self._array[0],
            self._array[1],
            self._array[2],
            self._array[3],
            self._array[4],
            self._array[5],
            self._array[6],
            self._array[7],
            self._array[8]
        )
    }
}

impl std::ops::Add<Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;
    fn add(self, other: Matrix3x3) -> Matrix3x3
    {
        let mut result = Matrix3x3::new();
        for i in 0..9 {
            result[i] = self[i] + other[i];
        }
        result
    }
}

impl std::ops::Add<&Matrix3x3> for &Matrix3x3 {
    type Output = Matrix3x3;
    fn add(self, other: &Matrix3x3) -> Matrix3x3 {
        *self + *other
    }
}

impl std::ops::Add<Matrix3x3> for &Matrix3x3 {
    type Output = Matrix3x3;
    fn add(self, other: Matrix3x3) -> Matrix3x3 {
        *self + other
    }
}

impl std::ops::Add<&Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;
    fn add(self, other: &Matrix3x3) -> Matrix3x3 {
        self + *other
    }
}

impl std::ops::Sub<Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;
    fn sub(self, other: Matrix3x3) -> Matrix3x3 {

        let mut result = Matrix3x3::new();
        for i in 0..9 {
            result[i] = self[i] - other[i];
        }
        return result
    }
}

impl std::ops::Sub<&Matrix3x3> for &Matrix3x3 {
    type Output = Matrix3x3;
    fn sub(self, other: &Matrix3x3) -> Matrix3x3 {
        *self - *other
    }
}

impl std::ops::Sub<Matrix3x3> for &Matrix3x3 {
    type Output = Matrix3x3;
    fn sub(self, other: Matrix3x3) -> Matrix3x3 {
        *self - other
    }
}

impl std::ops::Sub<&Matrix3x3> for Matrix3x3 {
    type Output = Matrix3x3;
    fn sub(self, other: &Matrix3x3) -> Matrix3x3 {
        self - *other
    }
}

impl Mul<f64> for Matrix3x3 {
    type Output = Matrix3x3;

    fn mul(mut self, scalar: f64) -> Matrix3x3 {
        self.mul_scalar(scalar);
        self
    }
}

impl Mul<Matrix3x3> for f64 {
    type Output = Matrix3x3;

    fn mul(self, mut matrix: Matrix3x3) -> Matrix3x3 {
        matrix.mul_scalar(self);
        matrix
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
    /// Implements the new() method for the Matrix3x3 struct, which creates and returns a new
    /// Matrix3x3 instance with all elements initialized to 0.
    ///
    /// # Arguments
    ///
    /// * self - the Matrix3x3 instance to create.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mat = Matrix3x3::new();
    ///
    /// assert_eq!(Matrix3x3::from_array([0.0, 0.0, 0.0,0.0, 0.0, 0.0,0.0, 0.0, 0.0,]), mat)
    /// ```
    pub fn new() -> Self {
        Self::default()
    }

    /// Constructs a new Matrix3x3 object from an array of 9 f64 elements.
    ///
    /// # Arguments
    ///
    /// * m - An array of 9 f64 elements, which are used to fill the matrix in column-major order.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let m = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    /// let mtx = Matrix3x3::from_array(m);
    ///
    /// assert_eq!(mtx[0], 1.0);
    /// assert_eq!(mtx[4], 5.0);
    /// assert_eq!(mtx[8], 9.0);
    /// ```
    pub fn from_array(m: [f64; 9]) -> Self {
        Matrix3x3 { _array: m }
    }


    //Add this matrix to another matrix.
    /// Add the values of the passed Matrix3x3 object rhs to the caller Matrix3x3 object self.
    ///
    /// # Arguments
    /// * rhs - A reference to the Matrix3x3 object whose values are to be added to the caller object.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// lhs.add(&rhs);
    /// let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    /// ```
    pub fn add(&mut self, rhs: &Matrix3x3) {
        let lhs = self;
        for i in 0..9 {
            lhs[i] = lhs[i] + rhs[i];
        }
    }

    /// Subtract the given rhs matrix from the caller matrix.
    ///
    /// This function subtracts each element of the rhs matrix from the corresponding
    /// element of the caller matrix and stores the result in the caller matrix itself.
    ///
    /// # Arguments
    /// * rhs: a reference to a Matrix3x3 object that is subtracted from the caller matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// lhs.sub(&rhs);
    /// let expected = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    ///```
    pub fn sub(&mut self, rhs: &Matrix3x3) {
        let lhs = self;
        for i in 0..9 {
            lhs[i] = lhs[i] - rhs[i];
        }
    }

    /// Multiplies each element of the caller matrix object by a scalar value.
    ///
    /// # Arguments
    ///
    /// * rhs - A scalar value to multiply the matrix with.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let rhs = 2.0;
    /// lhs.mul_scalar(rhs);
    /// let expected = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    ///```
    pub fn mul_scalar(&mut self, rhs: f64) {
        let lhs = self;
        for i in 0..9 {
            lhs[i] = lhs[i] * rhs;
        }
    }

    /// Divide all elements of the caller matrix object with a scalar value rhs.
    ///
    /// # Arguments
    ///
    /// * rhs: The scalar value to divide with.
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]);
    /// let rhs = 2.0;
    /// lhs.div_scalar(rhs);
    /// let expected = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    /// ```
    pub fn div_scalar(&mut self, rhs: f64) {
        let lhs = self;

        for i in 0..9 {
            lhs[i] = lhs[i] / rhs;
        }
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
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Vec3;
    ///
    /// let lhs: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let mut rhs = Vec3::new(2.0, 2.0, 2.0);
    /// lhs.mul_vec_mut(&mut rhs);
    /// assert_eq!(12.0, rhs[0]);
    /// assert_eq!(30.0, rhs[1]);
    /// assert_eq!(48.0, rhs[2]);
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
    /// use bioshell_numerical::Matrix3x3;
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
    /// # Arguments
    ///
    /// * self - a reference to the calling Matrix3x3 object.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mat: Matrix3x3 = Matrix3x3::from_array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
    /// let _det: f64 = mat.det();
    /// assert_eq!(0.0, _det);
    ///```
    pub fn det(&self) -> f64 {
        let lhs = self;
        let a = lhs[0] * (lhs[4] * lhs[8] - lhs[5] * lhs[7]);
        let b = lhs[1] * (lhs[3] * lhs[8] - lhs[5] * lhs[6]);
        let c = lhs[2] * (lhs[3] * lhs[7] - lhs[4] * lhs[6]);
        let returns = a - b + c;
        return returns;
    }

    /// Calculates the inverse of the matrix in-place.
    ///
    /// # Arguments
    ///
    /// * self - A mutable reference to the caller matrix object.
    ///
    /// # Panics
    ///
    /// This function may panic if the matrix is not invertible, i.e. has a determinant of zero.
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
    /// lhs.inverse();
    /// let expected = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    ///```
    ///
    /// # Errors
    ///
    /// This function returns a Result type that wraps the inverse if it exists. If the matrix is not invertible,
    /// i.e. has a determinant of zero, it returns an error containing a string that explains the cause of the error.
    /// Here is an example of how to handle the error:
    ///
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut mat = Matrix3x3::new();
    ///
    /// // calculate the inverse and handle the error
    /// if let Err(e) = mat.inverse() {
    /// println!("Error: {}", e);
    /// }
    ///
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

    /// Constructs a 3x3 matrix from its individual elements cx_x, cx_y, cx_z, cy_x, cy_y, cy_z, cz_x, cz_y, cz_z.
    ///
    /// # Arguments
    ///
    /// * cx_x - The element at row 0 and column 0 of the matrix.
    /// * cx_y - The element at row 0 and column 1 of the matrix.
    /// * cx_z - The element at row 0 and column 2 of the matrix.
    /// * cy_x - The element at row 1 and column 0 of the matrix.
    /// * cy_y - The element at row 1 and column 1 of the matrix.
    /// * cy_z - The element at row 1 and column 2 of the matrix.
    /// * cz_x - The element at row 2 and column 0 of the matrix.
    /// * cz_y - The element at row 2 and column 1 of the matrix.
    /// * cz_z - The element at row 2 and column 2 of the matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::vec3::Vec3;
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let mat: Matrix3x3 =
    ///     Matrix3x3::from_values(cx.x, cx.y, cx.z, cy.x, cy.y, cy.z, cz.x, cz.y, cz.z);
    /// assert_eq!(cx.x, mat[0]); assert_eq!(cx.y, mat[1]); assert_eq!(cx.z, mat[2]);
    /// assert_eq!(cy.x, mat[3]); assert_eq!(cy.y, mat[4]); assert_eq!(cy.z, mat[5]);
    /// assert_eq!(cz.x, mat[6]); assert_eq!(cz.y, mat[7]); assert_eq!(cz.z, mat[8]);
    ///```
    pub fn from_values(
        cx_x: f64,
        cx_y: f64,
        cx_z: f64,
        cy_x: f64,
        cy_y: f64,
        cy_z: f64,
        cz_x: f64,
        cz_y: f64,
        cz_z: f64,
    ) -> Self {
        Self::from_array([cx_x, cx_y, cx_z, cy_x, cy_y, cy_z, cz_x, cz_y, cz_z])
    }

    /// Construct a Matrix3x3 from three column vectors, representing the x, y, and z axes of the matrix.
    ///
    /// # Arguments
    ///
    /// * cx - A reference to a Vec3 object representing the x axis of the matrix.
    /// * cy - A reference to a Vec3 object representing the y axis of the matrix.
    /// * cz - A reference to a Vec3 object representing the z axis of the matrix.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Vec3;
    ///
    /// let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let mat: Matrix3x3 = Matrix3x3::from_column_vectors(&cx, &cy, &cx);
    /// assert_eq!(cx.x, mat[0]); assert_eq!(cx.y, mat[1]); assert_eq!(cx.z, mat[2]);
    /// assert_eq!(cy.x, mat[3]); assert_eq!(cy.y, mat[4]); assert_eq!(cy.z, mat[5]);
    /// assert_eq!(cz.x, mat[6]); assert_eq!(cz.y, mat[7]); assert_eq!(cz.z, mat[8]);
    ///```
    pub fn from_column_vectors(cx: &Vec3, cy: &Vec3, cz: &Vec3) -> Self {
        Self::from_array([cx.x, cx.y, cx.z, cy.x, cy.y, cy.z, cz.x, cz.y, cz.z])
    }

    /// This method creates a Matrix3x3 from three Vec3 objects representing the columns of the matrix.
    ///
    /// # Arguments
    ///
    /// * rx - A reference to a Vec3 object representing the first column of the matrix.
    /// * ry - A reference to a Vec3 object representing the second column of the matrix.
    /// * rz - A reference to a Vec3 object representing the third column of the matrix.
    ///
    /// # Example
    ///```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Vec3;
    ///
    /// let cx: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cy: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let cz: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    /// let mat: Matrix3x3 = Matrix3x3::from_row_vectors(&cx, &cy, &cx);
    /// assert_eq!(cx.x, mat[0]); assert_eq!(cy.x, mat[1]); assert_eq!(cz.x, mat[2]);
    /// assert_eq!(cx.y, mat[3]); assert_eq!(cy.y, mat[4]); assert_eq!(cz.y, mat[5]);
    /// assert_eq!(cx.z, mat[6]); assert_eq!(cy.z, mat[7]); assert_eq!(cz.z, mat[8]);
    ///```
    pub fn from_row_vectors(rx: &Vec3, ry: &Vec3, rz: &Vec3) -> Self {
        Self::from_array([rx.x, ry.x, rz.x, rx.y, ry.y, rz.y, rx.z, ry.z, rz.z])
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
    /// use bioshell_numerical::matrix::Matrix3x3;
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
    /// use bioshell_numerical::matrix::Matrix3x3;
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
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let lhs = Matrix3x3::identity();
    /// let rhs = Matrix3x3::from_values(1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0);
    /// assert_eq!(lhs, rhs);
    /// ```
    pub fn identity() -> Self
    {
        return Matrix3x3::from_values(  1.0, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0);
    }

    /// Computes the inverse of a 3x3 matrix using the Gauss-Jordan elimination method.
    ///
    /// # Arguments
    ///
    /// * mat - A 3x3 matrix whose inverse is to be computed.
    ///
    /// # Returns
    ///
    /// * The inverse of the input matrix as a Matrix3x3.
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    ///
    /// let mut lhs: Matrix3x3 =
    ///     Matrix3x3::from_array([0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0]);
    /// lhs.inverse();
    /// let expected = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    /// for (i, val) in expected.iter().enumerate() {
    ///     assert_eq!(lhs[i], *val);
    /// }
    ///```
    pub fn inverse_s(mat: Matrix3x3) -> Matrix3x3
    {
        let mut mat_copy = mat.clone();
        mat_copy.inverse();
        return mat_copy;
    }

    /// Multiplies a 3x3 matrix with a 3D vector and returns the resulting 3D vector.
    ///
    /// # Arguments
    ///
    /// * mat - A reference to a Matrix3x3 object.
    /// * vec - A reference to a Vec3 object.
    ///
    /// # Example
    /// ```rust
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Vec3;
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
    /// use bioshell_numerical::matrix::Matrix3x3;
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
        let rows = 3;
        let cols =3;
        for i in 0..rows {
            for j in 0..cols {
                result[i * cols + j] = a[i] * b[j];
            }
        }
        return result;
    }
}
