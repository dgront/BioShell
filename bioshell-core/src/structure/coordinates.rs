use std::ops::{Index, IndexMut, Range};

use bioshell_numerical::{Vec3};

/// Stateless immutable view of coordinates
pub struct CoordinatesView<'a> {  pub points: &'a  Coordinates, }


#[derive(Clone, Debug)]
pub struct Coordinates {
    box_len: f64,
    box_len_half: f64,
    v: Vec<Vec3>,
    chains: Vec<Range<usize>>
}

macro_rules! wrap_coordinate_to_box {
    ($val:expr, $L:expr, $coord:expr) => {
        $coord = $val;
        if $coord > $L { $coord = $coord - $L}
        else {
            if $coord < 0.0 { $coord = $L + $coord}
        }
    }
}

macro_rules! closest_image {
    ($c1:expr, $c2:expr, $L: expr,$L2: expr, $delta:expr) => {
        $delta = $c1 - $c2;
        if $delta > 0.0 {
            if $delta > $L2 {$delta -= $L}
        } else {
            if $delta < -$L2 {$delta += $L}
        }
    }
}

impl Coordinates {

    pub fn new(n: usize) -> Coordinates {
        let mut v = if n > 0 { Vec::with_capacity(n) } else { Vec::new() };

        if n > 0 {
            let zero = Vec3::from_float(0.0);
            v.resize(n, zero);
        }
        let chains: Vec<Range<usize>> = vec![0..n];
        let l: f64 = 100000.0;
        return Coordinates {box_len: l, box_len_half: l/2.0, v, chains};
    }

    #[inline(always)]
    pub fn box_len(&self) -> f64 { self.box_len }

    #[inline(always)]
    pub fn set_box_len(&mut self, new_box_len: f64) {
        self.box_len = new_box_len;
        self.box_len_half = new_box_len / 2.0;
    }

    /// Returns the number of chains in this system
    pub fn count_chains(&self) -> usize { return self.chains.len(); }

    /// Provides the (half-open) range of atoms that belong to a given chain.
    ///
    /// Per rust convention used in ``std::ops::Range`` struct, the returned
    /// ``start..end`` range contains all atoms indexed by ``start <= idx < end``
    pub fn chain_range(&self, idx: usize) -> &Range<usize> { &self.chains[idx] }

    pub fn distance_square(&self, i: usize, j: usize) -> f64 {

        let mut d = self.v[i].x - self.v[j].x;
        let mut d2 = d * d;
        d = self.v[i].y - self.v[j].y;
        d2 += d * d;
        d = self.v[i].z - self.v[j].z;
        d2 += d * d;
        return d2;
    }

    pub fn closest_distance_square(&self, i: usize, j: usize) -> f64 {

        let mut d:f64;
        closest_image!(self.v[i].x, self.v[j].x, self.box_len, self.box_len_half, d);
        let mut d2 = d * d;
        closest_image!(self.v[i].y, self.v[j].y, self.box_len, self.box_len_half, d);
        d2 += d * d;
        closest_image!(self.v[i].z, self.v[j].z, self.box_len, self.box_len_half, d);

        return d2 + d*d;
    }

    pub fn closest_distance_square_to_vec(&self, i: usize, v: &Vec3) -> f64 {

        let mut d:f64;
        closest_image!(self.v[i].x, v.x, self.box_len, self.box_len_half, d);
        let mut d2 = d * d;
        closest_image!(self.v[i].y, v.y, self.box_len, self.box_len_half, d);
        d2 += d * d;
        closest_image!(self.v[i].z, v.z, self.box_len, self.box_len_half, d);

        return d2 + d*d;
    }

    pub fn delta_x(&self, i: usize, x: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.v[i].x,x, self.box_len, self.box_len_half, d);
        d
    }

    pub fn delta_y(&self, i: usize, y: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.v[i].y, y, self.box_len, self.box_len_half, d);
        d
    }

    pub fn delta_z(&self, i: usize, z: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.v[i].z, z, self.box_len, self.box_len_half, d);
        d
    }

    pub fn size(&self) -> usize { return self.v.len(); }

    pub fn x(&self, i:usize) -> f64 { self.v[i].x }

    pub fn y(&self, i:usize) -> f64 { self.v[i].y }

    pub fn z(&self, i:usize) -> f64 { self.v[i].z }

    pub fn set_x(&mut self, i:usize, x: f64) {  wrap_coordinate_to_box!(x, self.box_len, self.v[i].x); }

    pub fn set_y(&mut self, i:usize, y: f64) {  wrap_coordinate_to_box!(y, self.box_len, self.v[i].y); }

    pub fn set_z(&mut self, i:usize, z: f64) {  wrap_coordinate_to_box!(z, self.box_len, self.v[i].z); }

    pub fn set(&mut self, i:usize, x: f64, y: f64, z: f64) {
        wrap_coordinate_to_box!(x, self.box_len, self.v[i].x);
        wrap_coordinate_to_box!(y, self.box_len, self.v[i].y);
        wrap_coordinate_to_box!(z, self.box_len, self.v[i].z);
    }

    pub fn add(&mut self, i:usize, x: f64, y: f64, z: f64) {
        wrap_coordinate_to_box!(self.v[i].x + x, self.box_len, self.v[i].x);
        wrap_coordinate_to_box!(self.v[i].y + y, self.box_len, self.v[i].y);
        wrap_coordinate_to_box!(self.v[i].z + z, self.box_len, self.v[i].z);
    }

    /// Copy coordinates of i-th atom from a given rhs coordinates
    /// This method (unlike set()) does not apply PBC. To the contrary, it assumes the two systems:
    ///  this and RHS have exactly the same simulation box geometry
    pub fn copy(&mut self, i:usize, rhs: &Coordinates) {
        self.v[i].x = rhs.v[i].x;
        self.v[i].y = rhs.v[i].y;
        self.v[i].z = rhs.v[i].z;
    }
}

impl Index<usize> for Coordinates {
    type Output = Vec3;
    fn index(&self, i: usize) -> &Vec3 {
        &self.v[i]
    }
}

impl IndexMut<usize> for Coordinates {
    fn index_mut(&mut self, i: usize) -> &mut Vec3 {
        &mut self.v[i]
    }
}
