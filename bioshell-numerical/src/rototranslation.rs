use crate::matrix::Matrix3x3;
use crate::vec3::Vec3;
use std::fmt;

//region struct Rototranslation
pub struct Rototranslation {
    _origin: Vec3,
    _cos_theta: f64,
    _sin_theta: f64,
    _axis: Vec3,
    _rotation_matrix: Matrix3x3,
    _inverse_rotation_matrix: Matrix3x3,
}
//endregion

//region impl Debug
impl fmt::Debug for Rototranslation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rototranslation")
            .field("rotation_matrix", &self._rotation_matrix)
            .finish()
    }
}
//endregion

impl Rototranslation {
    //region pub fn around_axis(origin: Vec3, start: Vec3, end: Vec3, angle_rad: f64) -> Rototranslation
    pub fn around_axis(origin: Vec3, start: Vec3, end: Vec3, angle_rad: f64) -> Rototranslation {
        let mut axis = end - start;
        axis.normalize();
        let cos_theta = angle_rad.cos();
        let sin_theta = angle_rad.sin();

        let u_dot = Vec3::outer_product_mat(axis, axis);
        let u_cross = Matrix3x3::from_values(
            0.0, -axis.z, axis.y,
            axis.z, 0.0, -axis.x,
            -axis.y, axis.x, 0.0
        );

        let u_identity = Matrix3x3::identity();

        let rotation_matrix = cos_theta * u_identity
            + sin_theta * u_cross
            + (1.0 - cos_theta) * u_dot;

        let inverse_rotation_matrix = cos_theta * u_identity.clone()
            - sin_theta * u_cross.clone()
            + (1.0 - cos_theta) * u_dot.clone();

        return Rototranslation {
            _origin:origin,
            _cos_theta:cos_theta,
            _sin_theta:sin_theta,
            _axis:axis,
            _rotation_matrix:rotation_matrix,
            _inverse_rotation_matrix:inverse_rotation_matrix,
        };
    }
    //endregion

    //region pub fn apply_mut(&self, vector: &Vec3) -> Vec3
    pub fn apply_mut(&self, vector: &Vec3) -> Vec3
    {
        let transformed_vector = *vector - self._origin;
        let temp = self._rotation_matrix;
        let transformed_vector = Vec3::transform(transformed_vector, temp);
        return transformed_vector + self._origin;
    }
    //endregion

    //region pub fn apply_inverse_mut(&self, v: &mut Vec3)
    pub fn apply_inverse_mut(&self, vector: &mut Vec3) -> Vec3  {
        let transformed_vector = *vector - self._origin;
        let temp = self._inverse_rotation_matrix;
        let transformed_vector = Vec3::transform(transformed_vector, temp);
        return transformed_vector + self._origin;
    }
    //endregion

    //region default functions
    pub fn apply_inverse(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_inverse_mut(&mut v);
        return v;
    }

    pub fn apply(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_mut(&mut v);
        return v;
    }
    //endregion
}
