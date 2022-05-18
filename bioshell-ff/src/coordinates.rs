use std::ops::{Index, IndexMut};

use bioshell_numerical::{Vec3};

#[derive(Clone)]
pub struct Coordinates {
    v: Vec<Vec3>,
}

impl Coordinates {
    pub fn new(n: usize) -> Coordinates {
        let mut v = Vec::with_capacity(n);
        let mut zero = Vec3::from_float(0.0);
        v.resize(n, zero);

        return Coordinates {v};
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
