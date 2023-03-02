use crate::Vec3;

#[derive(Clone)]
/// Linear transformations from one vector in space to another.
pub struct Rototranslation {
    /// Rotation matrix represented row-wise
    pub m: [f64; 9],
    /// translation vector
    pub v: Vec3
}

/// Default rototranslation does not change anything.
///
/// Rotation matrix is set to a unit matrix, translation vector set to zero.
impl Default for Rototranslation {
    fn default() -> Self {
        Rototranslation{ m: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], v: Vec3::new(0.0, 0.0, 0.0)}
    }
}

impl Rototranslation {
    /// Creates a rototranslation from given column vectors.
    ///
    /// This function takes care for the normalisation of the input vectors.
    pub fn from_column_vectors(cx: &Vec3, cy: &Vec3, cz: &Vec3, center: &Vec3) -> Rototranslation {
        Rototranslation::default()
    }

    /// Creates a rototranslation from given row vectors.
    ///
    /// This function takes care for the normalisation of the input vectors.
    pub fn from_row_vectors(rx: &Vec3, ry: &Vec3, rz: &Vec3, center: &Vec3) -> Rototranslation {
        Rototranslation::default()
    }

    /// Creates a transformation that rotate around an axis
    ///
    /// # Arguments
    /// * `center` - rotation center
    /// * `begin` - rotation axis starts here
    /// * `end` - second point necessary to define the axis of rotation
    /// * `angle` - angle of rotation
    pub fn around_axis(center: &Vec3, begin: &Vec3, end: &Vec3, angle: f64) -> Rototranslation {
        Rototranslation::default()
    }

    /// Returns a transformed copy of a given vector
    pub fn apply(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_mut(&mut v);
        return v;
    }

    /// Apply this rototranslation to a given vector
    pub fn apply_mut(&self, v: &mut Vec3) {

    }

    /// Apply the inverse of this rototranslation to a given vector.
    /// Returns a transformed copy of a given vector
    pub fn apply_inverse(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_inverse_mut(&mut v);
        return v;
    }

    /// Apply the inverse of this rototranslation to a given vector
    pub fn apply_inverse_mut(&self, v: &mut Vec3) {

    }
}