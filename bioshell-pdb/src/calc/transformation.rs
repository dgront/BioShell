use std::fmt;
use crate::calc::{Matrix3x3, Vec3};

/// Rotation-translation operation in 3D
pub struct Rototranslation {
    _origin: Vec3,
    translate_back: bool,
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
            translate_back: true,
            _rotation_matrix: u_rot,
            _inverse_rotation_matrix: inv,
        };
    }

    /// Rototranslation transforming to a local coordinate system of three atoms.
    ///
    /// The axes of the local coordinate system are defined by three points ``a``, ``b`` and ``c``,
    /// two versors are defined:
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
    /// let rot_m = rot.rotation_matrix();
    /// assert_delta!(rot_m.elem(0,0), 1.0, 0.000001);
    /// # assert_delta!(rot_m.elem(1,1), 1.0, 0.000001);
    /// # assert_delta!(rot_m.elem(2,2), 1.0, 0.000001);
    /// # println!("{}", rot);
    /// # let expected = r#"[  1.000  -0.000   0.000 ]     | vx -  2.000 |
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

        let mut n_to_ca = b.clone();  // --- a1 -> b vector
        n_to_ca -= a;
        n_to_ca.normalize();

        let mut ca_to_c = c.clone();  // --- a1 -> b vector
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

        let u_rot = Matrix3x3::from_column_vectors(&tx, &ty, &tz);
        let mut inv = u_rot.clone();

        inv.inverse();
        return Rototranslation {
            _origin: b.clone(),
            translate_back: false,
            _rotation_matrix: u_rot,
            _inverse_rotation_matrix: inv,
        };
    }

    /// Provides read-only access to the rotation matrix of this [`Rototranslation`](Rototranslation)
    pub fn rotation_matrix(&self) -> &Matrix3x3 { &self._rotation_matrix }

    /// Provides read-only access to translation vector of this [`Rototranslation`](Rototranslation)
    pub fn center(&self) -> &Vec3 { &self._origin }

    /// Sets the translation vector of this [`Rototranslation`](Rototranslation)
    pub fn set_center(&mut self, center: &Vec3) {
        self._origin.set(center);
    }

    /// Defines whether this [`Rototranslation`](Rototranslation) moves points back to the original center.
    ///
    /// By default, this [`Rototranslation`](Rototranslation) moves points to the center of the coordinate system,
    /// applies the rotation matrix and then moves points back to the original center.
    /// Use ```if_translate_back(false)`` to leave transformed points at the center of the coordinate system.
    pub fn if_translate_back(&mut self, flag: bool) {
        self.translate_back = flag
    }

    /// Transpose the rotation matrix of this [`Rototranslation`](Rototranslation)
    ///
    /// This effects the inverse rotation.
    pub fn transpose(&mut self) {
        self._rotation_matrix.transpose();
    }

    pub fn apply_mut(&self, vector: &mut Vec3) {
        *vector -= &self._origin;
        self._rotation_matrix.mul_vec_mut(vector);
        if self.translate_back { *vector += &self._origin; }
    }

    pub fn apply_inverse_mut(&self, vector: &mut Vec3) {
        if self.translate_back { *vector += &self._origin; }
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
        let cx = self._origin.x;
        let cy = self._origin.y;
        let cz = self._origin.z;

        if self.translate_back {
            print_row!(f, self._rotation_matrix, 0, " ", "vx", cx, format!("   | {:>6.3} |", cx))?;
            print_row!(f, self._rotation_matrix, 1, "*", "vy", cy, format!(" + | {:>6.3} |", cy))?;
            print_row!(f, self._rotation_matrix, 2, " ", "vz", cz, format!("   | {:>6.3} |", cz))?;
        } else {
            print_row!(f, self._rotation_matrix, 0, " ", "vx", cx, "")?;
            print_row!(f, self._rotation_matrix, 1, "*", "vy", cy, "")?;
            print_row!(f, self._rotation_matrix, 2, " ", "vz", cz, "")?;
        }

        Ok(())
    }
}