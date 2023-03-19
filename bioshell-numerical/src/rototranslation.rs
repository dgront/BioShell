use crate::vec3::Vec3;
use crate::matrix::Matrix3x3;

pub struct Rototranslation
{
    rotation_matrix: Matrix3x3,
    translation_vec: Vec3,
}

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

    pub fn show(&self)
    {
        self.rotation_matrix.show();
        self.translation_vec.show();
    }

    pub fn rotation_matrix(&self) -> &Matrix3x3
    {
        return &self.rotation_matrix;
    }
    pub fn center_vec(&self) -> &Vec3
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

    /// Creates a transformation that rotate around an axis
    ///
    /// # Arguments
    /// * `center` - rotation center
    /// * `begin` - rotation axis starts here
    /// * `end` - second point necessary to define the axis of rotation
    /// * `angle` - angle of rotation
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

    pub fn apply_mut(&self, v: &mut Vec3)
    {
        let temp_x = v.x * self.rotation_matrix[0] + v.y * self.rotation_matrix[1] + v.z * self.rotation_matrix[2];
        let temp_y = v.x * self.rotation_matrix[3] + v.y * self.rotation_matrix[4] + v.z * self.rotation_matrix[5];
        let temp_z = v.x * self.rotation_matrix[6] + v.y * self.rotation_matrix[7] + v.z * self.rotation_matrix[8];

        v.x = temp_x;
        v.y = temp_y;
        v.z = temp_z;
    }

    /// Returns a transformed copy of a given vector
    pub fn apply(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_mut(&mut v);
        return v;
    }

    pub fn apply_inverse_mut(&self, v: &mut Vec3)
    {
        let temp_x = v.x * self.rotation_matrix[0] + v.y * self.rotation_matrix[3] + v.z * self.rotation_matrix[6];
        let temp_y = v.x * self.rotation_matrix[1] + v.y * self.rotation_matrix[4] + v.z * self.rotation_matrix[7];
        let temp_z = v.x * self.rotation_matrix[2] + v.y * self.rotation_matrix[5] + v.z * self.rotation_matrix[8];

        v.x = temp_x;
        v.y = temp_y;
        v.z = temp_z;
    }

    /// Apply the inverse of this rototranslation to a given vector.
    /// Returns a transformed copy of a given vector
    pub fn apply_inverse(&self, v: &Vec3) -> Vec3
    {
        let mut v = v.clone();
        self.apply_inverse_mut(&mut v);
        return v;
    }
}