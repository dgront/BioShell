use rand::distributions::Distribution;
use rand::thread_rng;
use rand_distr::Normal;
use std::fmt;
use std::ops::{Index, IndexMut};

/// 3D vector used to represent an interaction center in a simulation
///
/// Vec3 struct contains coordinates of a point, chain index and staple atom typing data. Its implementation
/// provides basic vector-type calculations.
///
/// The example below tests a few basic properties of a unit cube of with 1.0:
/// ```
/// # use bioshell_numerical::{dihedral_angle4, planar_angle3, Vec3};
/// let cube_points = [[0f64, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
///     [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0]];
/// let cube_vec: Vec<Vec3>  = cube_points.iter().map(|p| Vec3::new(p[0],p[1],p[2])).collect();
///
/// let mut center_expected = Vec3::new(0.5, 0.5, 0.5);
/// let mut center = Vec3::new(0.0, 0.0, 0.0);
/// for v in &cube_vec { center.add(&v) }
/// center.div(8.0);
///
/// assert!(center_expected.distance_to(&center) < 0.0001);
/// center.sub(&center_expected);
/// assert!(center.length() < 0.00001);
///
/// assert!((planar_angle3(&cube_vec[3], &cube_vec[0], &cube_vec[1]).to_degrees() - 90.0).abs() < 0.0001);
/// assert!((dihedral_angle4(&cube_vec[3], &cube_vec[0], &cube_vec[1], &cube_vec[5]).to_degrees() - 90.0).abs() < 0.0001);
/// ```
#[derive(Clone, Copy)]
pub struct Vec3 {
    /// the ``x`` coordinate of this vector
    pub x: f64,
    /// the ``y`` coordinate of this vector
    pub y: f64,
    /// the ``z`` coordinate of this vector
    pub z: f64,
    /// residue type assigned to this vector
    pub res_type: u8,
    /// atom type assigned to this vector
    pub atom_type: u8,
    /// index of a chain this atom belongs to
    pub chain_id: u16,
}

//region Indexing facility implementation.
impl Index<usize> for Vec3 {
    type Output = f64;

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
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Index out of range for Vec3"),
        }
    }
}
//endregion

impl fmt::Debug for Vec3 {
    /// Prints nicely 3D coordinates of a vector
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{:.3} {:.3} {:.3}]", self.x, self.y, self.z)
    }
}

impl PartialEq for Vec3 {
    fn eq(&self, other: &Self) -> bool {
        return self.x == other.x && self.y == other.y && self.z == other.z;
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

macro_rules! float3_operation {
    ($lhs:expr, $OP:tt, $x:expr, $y:expr, $z:expr) => {
        $lhs.x $OP $x;
        $lhs.y $OP $y;
        $lhs.z $OP $z;
    };
}

impl Vec3 {
    pub fn show(&self) {
        println!("[{},{},{}]", self.x, self.y, self.z);
        println!(">{},{},{}<", self.atom_type, self.chain_id, self.res_type);
    }
    /// Creates a new vector from given coordinates.
    ///
    /// ``res_type``, ``atom_type`` and ``chain_id`` are by default set to ``0``
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 {
            x: x,
            y: y,
            z: z,
            res_type: 0,
            atom_type: 0,
            chain_id: 0,
        }
    }

    /// Creates a new vector with all coordinates equal to a given value.
    ///
    /// ``res_type``, ``atom_type`` and ``chain_id`` are by default set to ``0``
    pub fn from_float(value: f64) -> Vec3 {
        Vec3 {
            x: value,
            y: value,
            z: value,
            res_type: 0,
            atom_type: 0,
            chain_id: 0,
        }
    }

    /// Assigns new content to this vector
    pub fn set(&mut self, v: &Vec3) {
        vec_operation!(self,v,=);
    }

    /// Turns self into the opposite vector
    /// Sum of a vector and its opposite should be zero:
    /// ```
    /// # use bioshell_numerical::Vec3;
    /// let v1 = Vec3::new(1.0, 2.0, 3.0);
    /// let mut v2 = v1.clone();
    /// v2.opposite();
    /// v2.add(&v1);
    /// assert!(v2.length() < 0.00001);
    /// ```
    pub fn opposite(&mut self) {
        self.x = -self.x;
        self.y = -self.y;
        self.z = -self.z;
    }

    /// Adds a vector to this vector
    pub fn add(&mut self, v: &Vec3) {
        vec_operation!(self,v,+=);
    }
    pub fn add_s(v1: &Vec3, v2: &Vec3) -> Vec3 {
        return Vec3::new(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
    }

    /// Subtracts a vector from this vector
    pub fn sub(&mut self, v: &Vec3) {
        vec_operation!(self,v,-=);
    }
    pub fn sub_s(v1: &Vec3, v2: &Vec3) -> Vec3 {
        return Vec3::new(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
    }

    /// Subtracts a vector from this vector
    pub fn mul(&mut self, f: f64) {
        scalar_operation!(self,f,*=);
    }

    /// Subtracts a vector from this vector
    pub fn div(&mut self, f: f64) {
        scalar_operation!(self,f,/=);
    }

    /// Adds float x, y, x values to this vector
    pub fn add3(&mut self, x: f64, y: f64, z: f64) {
        float3_operation!(self, +=, x, y, z);
    }

    /// Subtracts float x, y, x values from this vector
    pub fn sub3(&mut self, x: f64, y: f64, z: f64) {
        float3_operation!(self, -=, x, y, z);
    }

    /// Returns the squared length of this 3D vector
    pub fn length_squared(&self) -> f64 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }

    /// Returns a normalized copy of this vector
    /// ```
    /// # use bioshell_numerical::Vec3;
    ///
    /// let v = Vec3::new(3.0, 2.0, 1.0).normalized();
    /// assert!((v.length() - 1.0).abs() < 0.00001);
    /// ```
    pub fn normalized(&self) -> Vec3 {
        let mut v = self.clone();
        v.div(self.length());
        return v;
    }

    /// Normalizes this vector
    pub fn normalize(&mut self) {
        self.div(self.length());
    }

    /// Returns a length of this vector
    pub fn length(&self) -> f64 {
        return self.length_squared().sqrt();
    }

    /// Calculate a dot product of two vectors
    ///
    /// ```
    /// # use bioshell_numerical::Vec3;
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
    /// # use bioshell_numerical::Vec3;
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
    /// # use bioshell_numerical::Vec3;
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
            z: a.x * b.y - a.y * b.x,
            res_type: 0,
            atom_type: 0,
            chain_id: a.chain_id,
        };
    }
}

/// Calculates a planar angle between two vectors in 3D
pub fn planar_angle2(a: &Vec3, b: &Vec3) -> f64 {
    let v = Vec3::dot(a, b) as f64;
    return (v / (a.length() as f64 * b.length() as f64)).acos();
}

/// Calculates a planar angle of the a-b-c triangle in 3D
pub fn planar_angle3(a: &Vec3, b: &Vec3, c: &Vec3) -> f64 {
    let mut v1: Vec3 = Vec3::clone(a);
    v1.sub(b);
    let mut v2: Vec3 = Vec3::clone(c);
    v2.sub(b);
    return planar_angle2(&v1, &v2);
}

/// Calculates a dihedral angle of the a-b-c-d points in 3D
pub fn dihedral_angle4(a: &Vec3, b: &Vec3, c: &Vec3, d: &Vec3) -> f64 {
    let mut b0 = b.clone(); // b0 = -1.0*(b - a)
    b0.sub(a);
    b0.mul(-1.0);
    let mut b1 = c.clone(); // b1 = c - b
    b1.sub(b);
    b1.normalize(); // normalize b1
    let mut b2 = d.clone(); // b2 = d - c
    b2.sub(c);

    let mut v = b1.clone(); // v is the projection of b0 onto plane perpendicular to b1
    v.mul(-Vec3::dot(&b0, &b1));
    v.add(&b0);

    let mut w = b1.clone(); // v is the projection of of b2 onto plane perpendicular to b1
    w.mul(-Vec3::dot(&b2, &b1));
    w.add(&b2);

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
        z,
        res_type: center.res_type,
        atom_type: center.atom_type,
        chain_id: center.chain_id,
    };
}
