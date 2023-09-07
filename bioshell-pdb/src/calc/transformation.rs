use std::fmt;
use crate::calc::{Matrix3x3, Vec3};

pub struct Rototranslation {
    _origin: Vec3,
    _rotation_matrix: Matrix3x3,
    _inverse_rotation_matrix: Matrix3x3,
}

impl fmt::Debug for Rototranslation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Rototranslation")
            .field("rotation_matrix", &self._rotation_matrix)
            .finish()
    }
}

impl Rototranslation {
    /// Creates a transformation that rotates 3D points around a given axis
    ///
    /// cos_theta * u_identity
    //             + sin_theta * u_cross
    //             + (1.0 - cos_theta) * u_dot;
    pub fn around_axis(start: &Vec3, end: &Vec3, angle_rad: f64) -> Rototranslation {

        let mut axis = end.clone();
        axis -= start;
        axis.normalize();
        let cos_theta = angle_rad.cos();
        let sin_theta = angle_rad.sin();

        let mut u_dot = Vec3::outer(&axis, &axis);
        let mut u_cross = Matrix3x3::from_array(
            [0.0, -axis.z, axis.y,
            axis.z, 0.0, -axis.x,
            -axis.y, axis.x, 0.0]
        );

        let mut u_rot = Matrix3x3::identity();
        u_rot *= cos_theta;
        u_cross *= sin_theta;
        u_dot *= 1.0 - cos_theta;
        u_rot += &u_cross;
        u_rot += &u_dot;

        let mut inv = u_rot.clone();
        inv.inverse();
        return Rototranslation {
            _origin: start.clone(),
            _rotation_matrix: u_rot,
            _inverse_rotation_matrix: inv,
        };
    }

    pub fn rotation_matrix(&self) -> &Matrix3x3 { &self._rotation_matrix }

    pub fn apply_mut(&self, vector: &mut Vec3) {
        *vector -= &self._origin;
        self._rotation_matrix.mul_vec_mut(vector);
        *vector += &self._origin;
    }

    pub fn apply_inverse_mut(&self, vector: &mut Vec3) {
        *vector -= &self._origin;
        self._inverse_rotation_matrix.mul_vec_mut(vector);
        *vector += &self._origin;
    }

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
}
