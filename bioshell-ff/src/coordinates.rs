use std::ops::{Index, IndexMut};

use bioshell_numerical::{Vec3};

#[derive(Clone)]
pub struct Coordinates {
    v: Vec<Vec3>,
}

impl Coordinates {
    pub fn new(n: usize) -> Coordinates {
        let mut v = Vec::with_capacity(n);
        let zero = Vec3::from_float(0.0);
        v.resize(n, zero);

        return Coordinates {v};
    }

    pub fn distance_square(&self, i: usize, j: usize) -> f32 {

        let mut d = self.v[i].x - self.v[j].x;
        let mut d2 = d * d;
        d = self.v[i].y - self.v[j].y;
        d2 += d * d;
        d = self.v[i].z - self.v[j].z;
        d2 += d * d;
        return d2;
    }

    pub fn size(&self) -> usize { return self.v.len(); }
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
