use std::fmt;
use crate::calc::{Matrix3x3, Vec3};

/// Rotation-translation operation in 3D
pub struct Rototranslation {
    _translation: Vec3,
    _rotation_matrix: Matrix3x3,
    _inverse_translation: Vec3,
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

    pub fn new(rot: Matrix3x3, t: Vec3) -> Rototranslation {

        // inverse affine translation:
        // x = R^{-1}(x' - t) = R^{-1}x' + (-R^{-1}t)
        let mut inv = rot.clone();
        inv.inverse();
        let mut inv_t = Matrix3x3::mul_vec_s(&inv, &t);
        inv_t *= -1.0;

        return Rototranslation {
            _rotation_matrix: rot,
            _translation: t,
            _inverse_rotation_matrix: inv,
            _inverse_translation: inv_t,
        };
    }

    /// Creates a transformation that rotates 3D points around a given axis
    ///
    /// The rotation matrix is computed using the [Rodrigues' rotation formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula):
    /// cos_theta * u_identity
    ///             + sin_theta * u_cross
    ///             + (1.0 - cos_theta) * u_dot;
    pub fn around_axis(start: &Vec3, end: &Vec3, angle_rad: f64) -> Rototranslation {
        let mut axis = end.clone();
        axis -= start;
        axis.normalize();

        let cos_theta = angle_rad.cos();
        let sin_theta = angle_rad.sin();

        let mut u_dot = Vec3::outer(&axis, &axis);

        let mut u_cross = Matrix3x3::from_array([
            0.0,     -axis.z,  axis.y,
            axis.z,   0.0,    -axis.x,
            -axis.y,   axis.x,  0.0,
        ]);

        let mut u_rot = Matrix3x3::identity();
        u_rot *= cos_theta;

        u_cross *= sin_theta;
        u_dot *= 1.0 - cos_theta;

        u_rot += &u_cross;
        u_rot += &u_dot;

        // t = p - R p
        let rotated_start =  Matrix3x3::mul_vec_s(&u_rot, start);
        let mut t = start.clone();
        t -= &rotated_start;

        return Rototranslation::new(u_rot, t);
    }

    /// Rototranslation transforming to a local coordinate system of three atoms.
    ///
    /// The axes of the local coordinate system are defined by three points ``a``, ``b`` and ``c``.
    /// The point ``b`` is the origin of the local coordinate system,
    /// the first axis points along the ``v_ac`` vector,
    /// the third axis is the bisector of the angle formed by the ``a``, ``b`` and ``c`` points,
    ///and the second axis is defined as a cross product of the first and the third axes.
    /// The axes are computed as below:
    ///
    /// ```math
    /// \begin{matrix}
    /// v_{ab} = ||b-a||\\
    /// v_{bc} = ||c-b||\\
    /// \end{matrix}
    /// ```
    /// based on them the three axes of the local coordinate system are defined as
    /// ```math
    /// \begin{align*}
    /// v_x & = ||v_{ab}+v_{bc}||\\
    /// v_y & = v_z \times v_x\\
    /// v_z & = ||v_{ab}-v_{bc}||\\
    /// \end{align*}
    /// ```
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{assert_delta};
    /// use bioshell_pdb::calc::{Rototranslation, Vec3};
    /// let a = Vec3::new(0.0, 0.0, 0.0);
    /// let b = Vec3::new(2.0, 0.0, 2.0);
    /// let c = Vec3::new(4.0, 0.0, 0.0);
    /// let rot = Rototranslation::by_three_atoms(&a, &b, &c);
    /// let rot_m = rot.rotation();
    /// assert_delta!(rot_m.elem(0,0), 1.0, 0.000001);
    /// # assert_delta!(rot_m.elem(1,1), 1.0, 0.000001);
    /// # assert_delta!(rot_m.elem(2,2), 1.0, 0.000001);
    /// # println!("{}", rot);
    /// # let expected = r#"[  1.000   0.000   0.000 ]     | vx -  2.000 |
    /// # [  0.000   1.000   0.000 ]  *  | vy -  0.000 |
    /// # [  0.000   0.000   1.000 ]     | vz -  2.000 |
    /// # "#;
    /// # assert_eq!(format!("{}", rot), expected);
    /// ```
    ///
    /// The output printed should look like as below:
    /// ```txt
    /// [  1.000   0.000   0.000 ]     | vx -  2.000 |
    /// [  0.000   1.000   0.000 ]  *  | vy -  0.000 |
    /// [  0.000   0.000   1.000 ]     | vz -  2.000 |
    ///
    /// ```
    pub fn by_three_atoms(a: &Vec3, b: &Vec3, c: &Vec3) -> Rototranslation {

        let mut n_to_ca = b.clone();  // --- a -> b vector
        n_to_ca -= a;
        n_to_ca.normalize();

        let mut ca_to_c = c.clone();  // --- b -> c vector
        ca_to_c -= b;
        ca_to_c.normalize();

        let mut tz = n_to_ca.clone();
        tz -= &ca_to_c;
        tz.normalize();

        let mut tx = n_to_ca.clone();
        tx += &ca_to_c;
        tx.normalize();

        let mut ty = Vec3::cross(&n_to_ca, &ca_to_c);
        ty.normalize();

        let u_rot = Matrix3x3::from_row_vectors(&tx, &ty, &tz);

        return Rototranslation::new(u_rot, b.clone());
    }

    /// Provides read-only access to the rotation matrix of this [`Rototranslation`](Rototranslation)
    pub fn rotation(&self) -> &Matrix3x3 { &self._rotation_matrix }

    /// Provides read-only access to translation vector of this [`Rototranslation`](Rototranslation)
    pub fn translation(&self) -> &Vec3 { &self._translation }

    /// Sets the new translation vector of this [`Rototranslation`](Rototranslation)
    pub fn set_translation(&mut self, t: &Vec3) {

        self._translation.set(t);
        // inverse affine translation: x = R^-1 x' - R^-1 t
        self._inverse_translation = Matrix3x3::mul_vec_s(&self._inverse_rotation_matrix, &self._translation);
        self._inverse_translation *= -1.0;
    }

    /// Applies the inverse of this transformation to the given vector, modifying it in place.
    pub fn apply_mut(&self, vector: &mut Vec3) {
        self._rotation_matrix.mul_vec_mut(vector);
        *vector += &self._translation;
    }

    /// Applies this transformation to the given vector, modifying it in place.
    pub fn apply_inverse_mut(&self, vector: &mut Vec3) {
        self._inverse_rotation_matrix.mul_vec_mut(vector);
        *vector += &self._inverse_translation;
    }

    /// Applies the inverse of this transformation to the given vector and returns the result as a new vector.
    ///
    /// The original vector is not modified.
    pub fn apply_inverse(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_inverse_mut(&mut v);
        return v;
    }

    /// Applies this transformation to the given vector and returns the result as a new vector.
    ///
    /// The original vector is not modified.
    pub fn apply(&self, v: &Vec3) -> Vec3 {
        let mut v = v.clone();
        self.apply_mut(&mut v);
        return v;
    }
    //
    // pub fn batch_apply(&self, v: &[Vec3]) -> Vec3 {
    //     let mut v = v.clone();
    //     self.apply_mut(&mut v);
    //     return v;
    // }
}


// Define a macro to simplify printing a single row
macro_rules! print_row {
    ($f:expr, $array:expr, $row:expr, $op:expr, $vx:expr, $tb:expr, $res:expr) => {
        writeln!(
            $f,
            "[ {:>6.3}  {:>6.3}  {:>6.3} ]  {}  | {} - {:>6.3} |{}",
            $array[$row * 3], $array[$row * 3 + 1], $array[$row * 3 + 2], $op, $vx, $tb, $res
        )
    };
}

/// Implement the Display trait for Rototranslation
impl fmt::Display for Rototranslation {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let cx = self._translation.x;
        let cy = self._translation.y;
        let cz = self._translation.z;

        print_row!(f, self._rotation_matrix, 0, " ", "vx", cx, "")?;
        print_row!(f, self._rotation_matrix, 1, "*", "vy", cy, "")?;
        print_row!(f, self._rotation_matrix, 2, " ", "vz", cz, "")?;

        Ok(())
    }
}