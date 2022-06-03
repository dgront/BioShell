use std::ops::{Index, IndexMut};

use std::io::stdout;
use std::io::Write;
use std::path::Path;
use std::fs::{File};

use bioshell_numerical::{Vec3};

#[derive(Clone)]
pub struct Coordinates {
    box_len: f32,
    box_len_half: f32,
    v: Vec<Vec3>,
}

impl Coordinates {
    pub fn new(n: usize) -> Coordinates {
        let mut v = Vec::with_capacity(n);
        let zero = Vec3::from_float(0.0);
        v.resize(n, zero);
        let l: f32 = 100000.0;
        return Coordinates {box_len: l, box_len_half: l/2.0, v};
    }

    #[inline(always)]
    pub fn box_len(&self) -> f32 { self.box_len }

    #[inline(always)]
    pub fn set_box_len(&mut self, new_box_len: f32) {
        self.box_len = new_box_len;
        self.box_len_half = new_box_len / 2.0;
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

    pub fn closest_distance_square(&self, i: usize, j: usize) -> f32 {

        let mut dx = self.v[i].x - self.v[j].x;
        if dx > self.box_len_half { dx -= self.box_len; }
        let mut d2 = dx * dx;
        dx = self.v[i].y - self.v[j].y;
        if dx > self.box_len_half { dx -= self.box_len; }
        d2 += dx * dx;
        dx = self.v[i].z - self.v[j].z;
        if dx > self.box_len_half { dx -= self.box_len; }
        d2 += dx * dx;

        return d2;
    }

    pub fn size(&self) -> usize { return self.v.len(); }

    pub fn x(&self, i:usize) -> f32 { self.v[i].x * self.box_len }

    pub fn y(&self, i:usize) -> f32 { self.v[i].x * self.box_len }

    pub fn z(&self, i:usize) -> f32 { self.v[i].x * self.box_len }
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

pub fn to_pdb(chain: &Coordinates, i_model: i16, out_fname: &str) {

    let mut out_writer = match out_fname {
        "" => Box::new(stdout()) as Box<dyn Write>,
        _ => {
            let path = Path::new(out_fname);
            Box::new(File::options().append(true).create(true).open(&path).unwrap()) as Box<dyn Write>
        }
    };

    out_writer.write(format!("MODEL    {i_model}\n").as_bytes()).ok();
    for i in 0..chain.size() {
        out_writer.write(format!("ATOM   {:4}{}  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C\n",
                                 i+1, " CA ", i+1, chain[i].x, chain[i].y, chain[i].z).as_bytes()).ok();
    }
    out_writer.write(b"ENDMDL\n").ok();
}