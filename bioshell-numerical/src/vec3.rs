use std::fmt;
use rand::distributions::Distribution;
use rand::thread_rng;
use rand_distr::Normal;

/// 3D vector used to represent interaction center in a simulation
#[derive(Clone)]
pub struct Vec3 {
    /// the ``x`` coordinate of this vector
    pub x: f64,
    /// the ``y`` coordinate of this vector
    pub y: f64,
    /// the ``z`` coordinate of this vector
    pub z: f64,
    /// residue type assigned to this vector
    pub res_type: i8,
    /// atom type assigned to this vector
    pub atom_type: i8,
    /// index of a chain this atom belongs to
    pub chain_id: i16
}

impl fmt::Debug for Vec3 {
    /// Prints nicely 3D coordinates of a vector
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{:.3} {:.3} {:.3}]", self.x, self.y, self.z)
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
    /// Creates a new vector from given coordinates.
    ///
    /// ``res_type``, ``atom_type`` and ``chain_id`` are by default set to ``0``
    pub fn new(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3 { x: x, y: y, z: z, res_type:0, atom_type: 0, chain_id: 0}
    }

    /// Creates a new vector with all coordinates equal to a given value.
    ///
    /// ``res_type``, ``atom_type`` and ``chain_id`` are by default set to ``0``
    pub fn from_float(value: f64) -> Vec3 {
        Vec3 {
            x: value,
            y: value,
            z: value,
            res_type:0, atom_type: 0, chain_id: 0
        }
    }

    /// Assigns new content to this vector
    pub fn set(&mut self, v: &Vec3) { vec_operation!(self,v,=); }

    /// Adds a vector to this vector
    pub fn add(&mut self, v: &Vec3) { vec_operation!(self,v,+=); }

    /// Subtracts a vector from this vector
    pub fn sub(&mut self, v: &Vec3) { vec_operation!(self,v,-=); }

    /// Subtracts a vector from this vector
    pub fn mul(&mut self, f: f64) { scalar_operation!(self,f,*=); }

    /// Subtracts a vector from this vector
    pub fn div(&mut self, f: f64) { scalar_operation!(self,f,/=); }

    /// Adds float x, y, x values to this vector
    pub fn add3(&mut self, x:f64, y:f64, z:f64) { float3_operation!(self, +=, x, y, z); }

    /// Subtracts float x, y, x values from this vector
    pub fn sub3(&mut self, x:f64, y:f64, z:f64) { float3_operation!(self, -=, x, y, z); }

    /// Returns the squared length of this 3D vector
    pub fn length_squared(&self) -> f64 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }

    /// Returns a normalized copy of this vector
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
    pub fn dot(a: &Vec3, b: &Vec3) -> f64 {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    /// Calculate vector product
    pub fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
        return Vec3 {
            x: a.y * b.z - a.z * b.y,
            y: a.z * b.x - a.x * b.z,
            z: a.x * b.y - a.y * b.x,
            res_type: 0,
            atom_type: 0,
            chain_id: a.chain_id
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
    return planar_angle2(&v1,&v2);
}


/// Calculates a dihedral angle of the a-b-c-d points in 3D
pub fn dihedral_angle4(a: &Vec3, b: &Vec3, c: &Vec3, d: &Vec3) -> f64 {

    let mut b0 = b.clone();             // b0 = -1.0*(b - a)
    b0.sub(a);
    b0.mul(-1.0);
    let mut b1 = c.clone();             // b1 = c - b
    b1.sub(b);
    b1.normalize();                            // normalize b1
    let mut b2 = d.clone();             // b2 = d - c
    b2.sub(c);

    let mut v = b1.clone();              // v is the projection of b0 onto plane perpendicular to b1
    v.mul(-Vec3::dot(&b0, &b1));
    v.add(&b0);

    let mut w = b1.clone();              // v is the projection of of b2 onto plane perpendicular to b1
    w.mul(-Vec3::dot(&b2, &b1));
    w.add(&b2);

    let x: f64 = Vec3::dot(&v, &w) as f64;
    let y: f64 = Vec3::dot(&Vec3::cross(&b1, &v), &w) as f64;

    return f64::atan2(y, x);
}

macro_rules! three_normal_rands {
    ($rng:expr) => {
        {
            let mut rng = thread_rng();
            let normal = Normal::new(0.0, 1.0).unwrap();
            (normal.sample(&mut rng), normal.sample(&mut rng), normal.sample(&mut rng))
        }
    };
}

/// Generates a random point on a unit sphere.
/// Coordinates of the newly generated point are returned as a tuple
pub fn random_unit_versor() -> (f64, f64, f64) {

    let (x, y, z) = three_normal_rands!(rand::thread_rng());
    let mut l: f64 = x * x + y * y + z * z;
    l = l.sqrt();
    return ((x/l) as f64, (y/l) as f64, (z/l) as f64);
}

/// Generates a random point that grows a chain in a random direction.
///
/// The newly generated point is randomly located on a sphere of a given `radius` and centered
/// on a given `center`
pub fn random_point_nearby(center: &Vec3, radius: f64) -> Vec3 {

    let (mut x, mut y, mut z) = three_normal_rands!(rand::thread_rng());
    let mut l: f64 = x * x + y * y + z * z;
    l = l.sqrt() / radius;
    x = x/l + center.x;
    y = x/l + center.y;
    z = x/l + center.z;
    return Vec3{x, y, z, res_type: center.res_type, atom_type: center.atom_type, chain_id: center.chain_id };
}