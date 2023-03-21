use std::fmt;
use crate::vec3::Vec3;
use crate::matrix::Matrix3x3;

/// Represents a rototranslation operation on a vector.
///
/// Internally the rotation-matrix elements are stored as a Matrix3x3 object,
/// and the translation vector is stored as a Vec3 object
///
/// # Example
/// ```rust
/// use bioshell_numerical::matrix::Matrix3x3;
/// use bioshell_numerical::Vec3;
///
/// let vx = Vec3::new(1.0, 0.0, 0.0);
///         let vy = Vec3::new(0.0, 1.0, 0.0);
///         let vz = Vec3::new(1.0, 0.0, 1.0);
///         let unit_vec = Vec3::new(1.0, 1.0, 1.0);
///         let another_vec = Vec3::new(10.0, 0.0, 10.0);
///         let unit_mtx = Matrix3x3::from_values(1.0, 0.0, 0.0,
///                                               0.0, 1.0, 0.0,
///                                               1.0, 0.0, 1.0);
///
///         let rototran = Rototranslation::new(unit_mtx, unit_vec);
///
///         rototran.apply(&another_vec);
///
///         assert_eq!(10.0, another_vec[0]);
///         assert_eq!(0.0, another_vec[1]);
///         assert_eq!(10.0, another_vec[2]);
/// ```
pub struct Rototranslation
{
    rotation_matrix: Matrix3x3,
    translation_vec: Vec3,
}

//region Implementation of Debug trait
impl fmt::Debug for Rototranslation
{
    /// Formats the struct for printing using the `debug_struct` macro and prints the values of `rotation_matrix` and `translation_vec`.
    ///
    /// # Arguments
    ///
    /// * `self` - A reference to the `Rototranslation` object being formatted
    /// * `f` - A mutable reference to a formatter, used to specify the format of the output.
    ///
    /// # Returns
    ///
    /// A `fmt::Result` indicating whether formatting succeeded or failed.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rototranslation")
            .field("rotation_matrix", &self.rotation_matrix)
            .field("translation_vec", &self.translation_vec)
            .finish()
    }
}
//endregion

impl Rototranslation
{
    pub fn new(m: Matrix3x3, v: Vec3) -> Self
    {
        Rototranslation
        {
            rotation_matrix: m,
            translation_vec: v
        }
    }

    //region properties/accessors
    pub fn rotation_matrix(&self) -> &Matrix3x3
    {
        return &self.rotation_matrix;
    }
    pub fn translation_vec(&self) -> &Vec3
    {
        return &self.translation_vec;
    }
    pub fn set_rotation_matrix(&mut self, mat: Matrix3x3)
    {
        self.rotation_matrix = mat;
    }
    pub fn set_rotation_center(&mut self, center: Vec3)
    {
        self.translation_vec = center;
    }
    //endregion

    //region pub fn around_axis(center: &Vec3, begin: &Vec3, end: &Vec3, angle: f64)
    /// Creates a transformation that rotates a vector around an axis defined by two points.
    ///
    /// This function returns a `Rototranslation` object that applies a rotation around an axis
    /// defined by two points (`begin` and `end`) and an angle (`angle`) with `center` as the center of rotation.
    ///
    /// # Arguments
    ///
    /// * `center` - The center point of the rotation.
    /// * `begin` - The starting point of the axis.
    /// * `end` - The ending point of the axis.
    /// * `angle` - The angle of rotation in radians.
    ///
    /// # Example
    ///
    /// ```rust
    /// let rotation_mat = matrix_lib::matrix::Matrix3x3::from_values(
    ///             1.0,0.0,0.0,
    ///             0.0,1.0,0.0,
    ///             0.0,0.0,1.0
    ///         );
    ///         let translation_vec = matrix_lib::vec3::Vec3::new(1.0,1.0,1.0);
    ///         let mut another_vec = Vec3::new(10.0, 1.0, 10.0);
    ///         let rot = Rototranslation::new(rotation_mat, translation_vec);
    ///
    ///         let center = Vec3::new(0.0, 0.0, 0.0);
    ///         let begin = Vec3::new(0.0, 0.0, 1.0);
    ///         let end = Vec3::new(0.0, 1.0, 0.0);
    ///         let angle = std::f64::consts::PI / 2.0;
    ///         let rot = Rototranslation::around_axis(&center, &begin, &end, angle);
    ///
    ///         assert_eq!(10.0, another_vec.x);
    ///         assert_eq!(1.0, another_vec.y);
    ///         assert_eq!(10.0, another_vec.z);
    /// ```
    pub fn around_axis(center: &Vec3, begin: &Vec3, end: &Vec3, angle: f64) -> Rototranslation
    {
        let cosa: f64 = angle.cos();
        let sina: f64 = angle.sin();
        let axis = Vec3::sub_s(end, begin);
        let x: f64 = axis.x;
        let y: f64 = axis.y;
        let z: f64 = axis.z;
        let x2: f64 = axis.x * axis.x;
        let y2: f64 = axis.y * axis.y;
        let z2: f64 = axis.z * axis.z;
        let vt: f64 = 1.0 - cosa;

        let rot_x_x = x2 * vt + cosa;
        let rot_x_y = x * y * vt - z * sina;
        let rot_x_z = x * z * vt + y * sina;

        let rot_y_x = x * y * vt + z * sina;
        let rot_y_y = y2 * vt + cosa;
        let rot_y_z = y * z * vt - x * sina;

        let rot_z_x = x * z * vt - y * sina;
        let rot_z_y = y * z * vt + x * sina;
        let rot_z_z = z2 * vt + cosa;

        let mat: Matrix3x3 = Matrix3x3::from_values(rot_x_x, rot_x_y, rot_x_z,
                                                    rot_y_x, rot_y_y, rot_y_z,
                                                    rot_z_x, rot_z_y, rot_z_z);
        return Self::new(mat, Vec3::new(center.x, center.y, center.z));
    }
    //endregion

    //region pub fn apply_mut(&self, v: &mut Vec3)
    /// Applies the rototranslation operation to a mutable vector in place.
    ///
    /// This method applies the rotation matrix and the translation vector to the input vector
    /// and stores the result back in the input vector. The operation is performed in place,
    /// meaning the input vector is modified.
    ///
    /// # Arguments
    ///
    /// * v - A mutable reference to a Vec3 object representing the vector to be transformed.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bioshell_numerical::Vec3;
    /// use bioshell_numerical::Rototranslation;
    ///
    /// let rotation_mat = matrix_lib::matrix::Matrix3x3::from_values(
    ///             1.0,0.0,0.0,
    ///             0.0,1.0,0.0,
    ///             0.0,0.0,1.0
    ///         );
    ///         let translation_vec = matrix_lib::vec3::Vec3::new(1.0,1.0,1.0);
    ///         let mut another_vec = Vec3::new(10.0, 1.0, 10.0);
    ///
    ///         let rot = Rototranslation::new(rotation_mat, translation_vec);
    ///
    ///         rot.apply_mut(&mut another_vec);
    ///
    ///         assert_eq!(10.0, another_vec.x);
    ///         assert_eq!(1.0, another_vec.y);
    ///         assert_eq!(10.0, another_vec.z);
    /// ```
    pub fn apply_mut(&self, v: &mut Vec3)
    {
        let temp_x = v.x * self.rotation_matrix[0] + v.y * self.rotation_matrix[1] + v.z * self.rotation_matrix[2];
        let temp_y = v.x * self.rotation_matrix[3] + v.y * self.rotation_matrix[4] + v.z * self.rotation_matrix[5];
        let temp_z = v.x * self.rotation_matrix[6] + v.y * self.rotation_matrix[7] + v.z * self.rotation_matrix[8];

        v.x = temp_x;
        v.y = temp_y;
        v.z = temp_z;
    }
    //endregion

    //region pub fn apply(&self, v: &Vec3) -> Vec3
    /// Applies the rototranslation to a vector.
    ///
    /// # Arguments
    ///
    /// * v - A reference to the vector to apply the rototranslation to
    ///
    /// # Example
    ///
    /// ```rust
    /// use bioshell_numerical::Vec3;
    /// use bioshell_numerical::Rototranslation;
    ///
    /// let unit_matrix =
    ///             Matrix3x3::from_values(1.0, 0.0, 0.0,
    ///                                     0.0, 1.0, 0.0,
    ///                                     0.0, 0.0, 1.0);
    ///         let vector = Vec3::new(10.0, 0.0, 10.0);
    ///
    ///         let rot = Rototranslation::new(unit_matrix, vector);
    ///         let rotated_vector = rot.apply(&vector);
    ///
    ///         assert_eq!(rotated_vector.x, 10.0);
    ///         assert_eq!(rotated_vector.y, 0.0);
    ///         assert_eq!(rotated_vector.z, 10.0);
    /// ```
    pub fn apply(&self, v: &Vec3) -> Vec3
    {
        let mut v = v.clone();
        self.apply_mut(&mut v);
        return v;
    }
    //endregion

    //region pub fn apply_inverse_mut(&self, v: &mut Vec3)
    /// Applies the inverse of the Rototranslation object to a mutable reference to a Vec3 object.
    ///
    /// This function applies the inverse transformation to the input vector in place. It first translates
    /// the vector by the negation of the translation vector and then rotates it by the inverse of the
    /// rotation matrix.
    ///
    /// # Arguments
    ///
    /// * v - a mutable reference to a Vec3 object that will be transformed in place
    ///
    /// # Example
    ///
    /// ```rust
    /// use bioshell_numerical::Vec3;
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Rototranslation;
    ///
    /// let unit_matrix = Matrix3x3::from_values(  1.0, 0.0, 0.0,
    ///                                                             0.0, 1.0, 0.0,
    ///                                                             0.0, 0.0, 1.0);
    ///         let trans_vec = Vec3::new(1.0, 1.0, 1.0);
    ///         let rt = Rototranslation::new(unit_matrix, trans_vec);
    ///         let mut another_vec = Vec3::new(1.0, 2.0, 3.0);
    ///
    ///         rt.apply_inverse_mut(&mut another_vec);
    ///
    ///         assert_eq!(Vec3::new(1.0, 2.0, 3.0), another_vec);
    ///
    pub fn apply_inverse_mut(&self, v: &mut Vec3)
    {
        let temp_x = v.x * self.rotation_matrix[0] + v.y * self.rotation_matrix[3] + v.z * self.rotation_matrix[6];
        let temp_y = v.x * self.rotation_matrix[1] + v.y * self.rotation_matrix[4] + v.z * self.rotation_matrix[7];
        let temp_z = v.x * self.rotation_matrix[2] + v.y * self.rotation_matrix[5] + v.z * self.rotation_matrix[8];

        v.x = temp_x;
        v.y = temp_y;
        v.z = temp_z;
    }
    //endregion

    //region pub fn apply_inverse(&self, v: &Vec3) -> Vec3
    /// Applies the inverse of the rototranslation to a vector in place.
    ///
    /// The inverse of the rototranslation is obtained by inverting the rotation matrix and
    /// by reversing the translation direction.
    ///
    /// # Arguments
    ///
    /// * v - The vector to which the inverse of the rototranslation is applied.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bioshell_numerical::Vec3;
    /// use bioshell_numerical::matrix::Matrix3x3;
    /// use bioshell_numerical::Rototranslation;
    ///
    /// let unit_matrix = Matrix3x3::from_values(  1.0, 0.0, 0.0,
    ///                                                    0.0, 1.0, 0.0,
    ///                                                    0.0, 0.0, 1.0);
    ///         let trans_vec = Vec3::new(1.0, 1.0, 1.0);
    ///         let rt = Rototranslation::new(unit_matrix, trans_vec);
    ///         let another_vec = Vec3::new(1.0, 2.0, 3.0);
    ///
    ///         let vec_out = rt.apply_inverse(&another_vec);
    ///
    ///         assert_eq!(Vec3::new(1.0, 2.0, 3.0), vec_out);
    /// ```
    pub fn apply_inverse(&self, v: &Vec3) -> Vec3
    {
        let mut v = v.clone();
        self.apply_inverse_mut(&mut v);
        return v;
    }
    //endregion
}