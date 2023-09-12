use rand::distributions::Distribution;
use rand::thread_rng;
use rand_distr::Normal;
use std::fmt;
use std::fmt::{Display, Formatter};
use std::ops::{Index, IndexMut, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::calc::Matrix3x3;

/// 3D vector used to manipulate with atomic coordinates.
///
/// [`Vec3`] struct contains 3D coordinates; it is used to store the location of a [`PdbAtom`](crate::PdbAtom).
/// The [`Vec3`] struct implements also a few operators, such as `-=` or `+=` to facilitate vector arithmetics.
///
/// The example below tests a few basic properties of a unit cube of with 1.0, which looks as below:
/// ```text
///   h----g             The vertex a is placed at [0, 0, 0]
///  /    /|
/// e----f |             The vertex d is placed at [0, 1, 0] (hidden behind the front face abfe)
/// |    | c
/// |    |/
/// a----b
/// ```
/// ```
/// # use bioshell_pdb::calc::{dihedral_angle4, planar_angle3, Vec3};
/// let cube_points = [[0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]];
/// let cube_vec: Vec<Vec3>  = cube_points.iter().map(|p| Vec3::new(p[0],p[1],p[2])).collect();
/// let mut center_expected = Vec3::new(0.5, 0.5, 0.5);
/// let mut center = Vec3::from_float(0.0);
/// for v in &cube_vec { center += &v }
/// center /= 8.0;
///
/// assert!(center_expected.distance_to(&center) < 0.0001);
/// center -= &center_expected;
/// assert!(center.length() < 0.00001);
///
/// // unpack the vertices to allow calling them by letters as on the diagram above
/// let [a, b, c, d, e, f, g, h] = <[Vec3; 8]>::try_from(cube_vec).ok().unwrap();
/// assert!((planar_angle3(&a, &b, &f).to_degrees() - 90.0).abs() < 0.0001);
/// assert!((dihedral_angle4(&e, &a, &b, &c).to_degrees() + 90.0).abs() < 0.0001);
/// ```
#[derive(Clone, Copy, Default)]
pub struct Vec3 {
    /// the ``x`` coordinate of this vector
    pub x: f64,
    /// the ``y`` coordinate of this vector
    pub y: f64,
    /// the ``z`` coordinate of this vector
    pub z: f64
}

impl Index<usize> for Vec3 {
    type Output = f64;

    /// Indexing operator provides access to X, Y, Z components of a vector
    fn index(&self, index: usize) -> &f64 {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Index out of range for Vec3"),
        }
    }
}

impl IndexMut<usize> for Vec3 {
    /// Indexing operator provides mutable access to X, Y, Z components of a vector
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Index out of range for Vec3"),
        }
    }
}

macro_rules! vec_operation {
    ($lhs:expr, $rhs:expr, $OP:tt ) => {
        $lhs.x $OP $rhs.x;
        $lhs.y $OP $rhs.y;
        $lhs.z $OP $rhs.z;
    };
}

macro_rules! scalar_operation {
    ($lhs:expr, $rhs:expr, $OP:tt ) => {
        $lhs.x $OP $rhs;
        $lhs.y $OP $rhs;
        $lhs.z $OP $rhs;
    };
}

impl SubAssign<&Vec3> for Vec3 {
    /// Performs the `-=` operation.
    ///
    /// # Example
    ///
    /// ```
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v0 = Vec3::new(1.0, 2.0, 3.0);
    /// let mut v1 = Vec3::new(1.0, 2.0, 3.0);
    /// v0 -= &v1;
    /// assert!((v0.x).abs() < 0.000001);
    /// ```
    fn sub_assign(&mut self, other: &Vec3) {
        vec_operation!(self, other, -=);
    }
}

impl AddAssign<&Vec3> for Vec3 {
    /// Performs the `+=` operation.
    ///
    /// # Example
    ///
    /// ```
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v0 = Vec3::new(1.0, 2.0, 3.0);
    /// let mut v1 = Vec3::new(1.0, 1.0, 1.0);
    /// v0 += &v1;
    /// assert!((v0.x - 2.0).abs() < 0.000001);
    /// ```
    fn add_assign(&mut self, other: &Vec3) {
        vec_operation!(self, other, +=);
    }
}

impl MulAssign<f64> for Vec3 {
    /// Performs the `*=` operation that multiplies this vector by a constant.
    ///
    /// ```
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v0 = Vec3::new(1.0, 2.0, 3.0);
    /// v0 *= 2.0;
    /// assert!((v0.x - 2.0).abs() < 0.000001);
    /// ```
    fn mul_assign(&mut self, rhs: f64) {
        scalar_operation!(self, rhs, *=);
    }
}

impl DivAssign<f64> for Vec3 {
    /// Performs the `/=` operation that divides this vector by a constant.
    ///
    /// ```
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v0 = Vec3::new(1.0, 2.0, 4.0);
    /// v0 /= 2.0;
    /// assert!((v0.z - 2.0).abs() < 0.000001);
    /// ```
    fn div_assign(&mut self, rhs: f64) {
        scalar_operation!(self, rhs, /=);
    }
}

impl fmt::Debug for Vec3 {
    /// Debug formatting of a Vec3 prints all its fields, e.g.
    /// ```rust
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v = Vec3::new(0.0, 1.0, 2.0);
    /// assert_eq!(format!("{:?}",v), "[0.000 1.000 2.000]");
    /// ```
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "[{:.3} {:.3} {:.3}]", self.x, self.y, self.z)
    }
}

impl Display for Vec3 {
    /// Prints X Y Z coordinates of a given 3D vector
    /// ```rust
    /// use bioshell_pdb::calc::Vec3;
    /// let mut v = Vec3::new(0.0, 1.0, 2.0);
    /// assert_eq!(format!("{}",v), "0.000 1.000 2.000");
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:.3} {:.3} {:.3}", self.x, self.y, self.z)
    }
}

impl PartialEq for Vec3 {
    /// Two vectors are equal if the values of all three coordinates are equal
    fn eq(&self, other: &Self) -> bool {
        return self.x == other.x && self.y == other.y && self.z == other.z;
    }
}

impl Vec3 {

    /// Creates a new vector from given coordinates.
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let v1 = Vec3::new(1.0, 2.0, 2.0);
    /// # assert_eq!(v1.length(), 3.0);
    /// ```
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x, y, z }
    }

    /// Creates a new vector with all coordinates equal to a given value.
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let zero_vec = Vec3::from_float(0.0);
    /// assert_eq!(zero_vec.length(), 0.0);
    /// ```
    pub fn from_float(value: f64) -> Vec3 {
        Vec3 { x: value, y: value, z: value }
    }

    /// Creates a new vector with given values
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let zero_vec = Vec3::from_array(&[2.0, 1.0, 2.0]);
    /// assert_eq!(zero_vec.length(), 3.0);
    /// ```
    pub fn from_array(values: &[f64; 3]) -> Vec3 {
        Vec3 { x: values[0], y: values[1], z: values[2] }
    }

    /// Assigns new content to this vector.
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let mut v1 = Vec3::new(1.0, 2.0, 2.0);
    /// v1.set(&Vec3::from_float(0.0));
    /// assert_eq!(v1.length(), 0.0);
    /// ```
    pub fn set(&mut self, v: &Vec3) {
        vec_operation!(self,v,=);
    }

    /// Add two vectors to creates a new one
    pub fn add_s(a: &Vec3, b: &Vec3) -> Vec3 {
        Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
    }

    /// Subtract two vectors to creates a new one
    pub fn sub_s(a: &Vec3, b: &Vec3) -> Vec3 {
        Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
    }

    /// Assigns new content to this vector.
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let mut v1 = Vec3::new(1.0, 2.0, 2.0);
    /// v1.set3(3.0, 0.0, 1.0);
    /// assert_eq!(v1.length(), 3.0);
    /// ```
    pub fn set3(&mut self, x: f64, y: f64, z: f64) {
        self.x = x;
        self.y = y;
        self.z = z;
    }

    /// Turns self into the opposite vector.
    ///
    /// Sum of a vector and its opposite should be zero:
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// let v1 = Vec3::new(1.0, 2.0, 3.0);
    /// let mut v2 = v1.clone();
    /// v2.opposite();
    /// v2 += &v1;
    /// assert!(v2.length() < 0.00001);
    /// ```
    pub fn opposite(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
    }

    /// Returns a length of this vector
    pub fn length(&self) -> f64 {
        return self.length_squared().sqrt();
    }

    /// Returns the squared length of this 3D vector
    pub fn length_squared(&self) -> f64 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }

    /// Returns a normalized copy of this vector
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    ///
    /// let v = Vec3::new(3.0, 2.0, 1.0).normalized();
    /// assert!((v.length() - 1.0).abs() < 0.00001);
    /// ```
    pub fn normalized(&self) -> Vec3 {
        let mut v = self.clone();
        v /= self.length();
        return v;
    }

    /// Normalizes this vector
    pub fn normalize(&mut self) { *self /= self.length(); }

    /// Calculate a dot product of two vectors
    ///
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// // let's try two ortogonal vectors
    /// let v1 = Vec3::new(3.0, 2.0, 1.0);
    /// let v2 = Vec3::new(-2.0, 3.0, 0.0);
    /// assert!((Vec3::dot(&v1, &v2)).abs() < 0.00001);
    /// ```
    pub fn dot(a: &Vec3, b: &Vec3) -> f64 {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /// Calculate the squared distance to another point
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// // Classic Pytagoras triangle with edges 3, 4 and 5
    /// let d = Vec3::new(3.0, 0.0, 0.0).distance_square_to(&Vec3::new(0.0, 4.0, 0.0));
    /// assert!((d-25.0).abs() < 0.00001);
    /// ```
    pub fn distance_square_to(&self, p: &Vec3) -> f64 {
        let mut d = self.x - p.x;
        let mut d2 = d * d;
        d = self.y - p.y;
        d2 += d * d;
        d = self.z - p.z;
        d2 += d * d;
        return d2;
    }
    /// Calculate the distance to another point
    pub fn distance_to(&self, p: &Vec3) -> f64 {
        self.distance_square_to(p).sqrt()
    }

    /// Calculate vector product
    /// ```
    /// # use bioshell_pdb::calc::Vec3;
    /// // multiply X and Y versors to get Z
    /// let x = Vec3::new(1.0, 0.0, 0.0);
    /// let y = Vec3::new(0.0, 1.0, 0.0);
    /// let z = Vec3::cross(&x, &y);
    /// assert!((z.z - 1.0).abs() < 0.0001);
    /// assert!((z.length() - 1.0).abs() < 0.0001);
    /// ```
    pub fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
        return Vec3 {
            x: a.y * b.z - a.z * b.y,
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x
        };
    }

    /// Calculate outer product of two vectors
    pub fn outer(lhs: &Vec3, rhs: &Vec3) -> Matrix3x3 {
        let mut data = [0.0; 9];
        data[0] = lhs.x * rhs.x;
        data[1] = lhs.x * rhs.y;
        data[2] = lhs.x * rhs.z;

        data[3] = lhs.y * rhs.x;
        data[4] = lhs.y * rhs.y;
        data[5] = lhs.y * rhs.z;

        data[6] = lhs.z * rhs.x;
        data[7] = lhs.z * rhs.y;
        data[8] = lhs.z * rhs.z;

        Matrix3x3::from_array(data)
    }
}

/// Calculates a planar angle between two vectors in 3D
pub fn planar_angle2(a: &Vec3, b: &Vec3) -> f64 {
    let v = Vec3::dot(a, b) as f64;
    return (v / (a.length() as f64 * b.length() as f64)).acos();
}

/// Calculates a planar angle of the a-b-c triangle in 3D
/// ```
/// use bioshell_pdb::calc::{planar_angle3, Vec3};
/// let va = Vec3::new(1.0, 0.0, 0.0);
/// let vb = Vec3::from_float(0.0);
/// let vc = Vec3::new(0.0, 1.0, 0.0);
/// assert!((planar_angle3(&va, &vb, &vc).to_degrees()-90.0)<0.001);
/// ```
pub fn planar_angle3(a: &Vec3, b: &Vec3, c: &Vec3) -> f64 {

    let mut v1: Vec3 = Vec3::clone(a);
    v1 -= b;
    let mut v2: Vec3 = Vec3::clone(c);
    v2 -= b;
    return planar_angle2(&v1, &v2);
}

/// Calculates a dihedral angle defined by the four a-b-c-d points in 3D
pub fn dihedral_angle4(a: &Vec3, b: &Vec3, c: &Vec3, d: &Vec3) -> f64 {

    let mut b0 = b.clone(); // b0 = -(b - a)
    b0 -= a;
    b0.opposite();
    let mut b1 = c.clone(); // b1 = c - b
    b1 -= b;
    b1.normalize(); // normalize b1
    let mut b2 = d.clone(); // b2 = d - c
    b2 -= c;

    let mut v = b1.clone(); // v is the projection of b0 onto plane perpendicular to b1
    v *= -Vec3::dot(&b0, &b1);
    v += &b0;

    let mut w = b1.clone(); // v is the projection of of b2 onto plane perpendicular to b1
    w *= -Vec3::dot(&b2, &b1);
    w += &b2;

    let x: f64 = Vec3::dot(&v, &w) as f64;
    let y: f64 = Vec3::dot(&Vec3::cross(&b1, &v), &w) as f64;

    return f64::atan2(y, x);
}

macro_rules! three_normal_rands {
    ($rng:expr) => {{
        let mut rng = thread_rng();
        let normal = Normal::new(0.0, 1.0).unwrap();
        (
            normal.sample(&mut rng),
            normal.sample(&mut rng),
            normal.sample(&mut rng),
        )
    }};
}

/// Generates a random point on a unit sphere.
///
/// Coordinates of the newly generated point are returned as a tuple
pub fn random_unit_versor() -> (f64, f64, f64) {
    let (x, y, z) = three_normal_rands!(rand::thread_rng());
    let mut l: f64 = x * x + y * y + z * z;
    l = l.sqrt();
    return ((x / l) as f64, (y / l) as f64, (z / l) as f64);
}

/// Generates a random point nearby a given location.
///
/// The newly generated point is randomly located on a sphere of a given `radius` and centered
/// on a given `center`
pub fn random_point_nearby(center: &Vec3, radius: f64) -> Vec3 {
    let (mut x, mut y, mut z) = three_normal_rands!(rand::thread_rng());
    let mut l: f64 = x * x + y * y + z * z;
    l = l.sqrt() / radius;
    x = x / l + center.x;
    y = y / l + center.y;
    z = z / l + center.z;
    return Vec3 {
        x,
        y,
        z
    };
}
