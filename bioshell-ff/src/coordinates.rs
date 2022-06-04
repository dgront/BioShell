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

        let mut d:f32;
        closest_image!(self.v[i].x, self.v[j].x, self.box_len, self.box_len_half, d);
        let mut d2 = d * d;
        closest_image!(self.v[i].y, self.v[j].y, self.box_len, self.box_len_half, d);
        d2 += d * d;
        closest_image!(self.v[i].z, self.v[j].z, self.box_len, self.box_len_half, d);

        return d2 + d*d;
    }

    pub fn delta_x(&self, i: usize, x: f32) -> f32 {
        let mut d: f32;
        closest_image!(self.v[i].x,x, self.box_len, self.box_len_half, d);
        d
    }

    pub fn delta_y(&self, i: usize, y: f32) -> f32 {
        let mut d: f32;
        closest_image!(self.v[i].y, y, self.box_len, self.box_len_half, d);
        d
    }

    pub fn delta_z(&self, i: usize, z: f32) -> f32 {
        let mut d: f32;
        closest_image!(self.v[i].z, z, self.box_len, self.box_len_half, d);
        d
    }

    pub fn size(&self) -> usize { return self.v.len(); }

    pub fn x(&self, i:usize) -> f32 { self.v[i].x }

    pub fn y(&self, i:usize) -> f32 { self.v[i].y }

    pub fn z(&self, i:usize) -> f32 { self.v[i].z }

    pub fn set_x(&mut self, i:usize, x: f32) {  wrap_coordinate_to_box!(x, self.box_len, self.v[i].x); }

    pub fn set_y(&mut self, i:usize, y: f32) {  wrap_coordinate_to_box!(y, self.box_len, self.v[i].y); }

    pub fn set_z(&mut self, i:usize, z: f32) {  wrap_coordinate_to_box!(z, self.box_len, self.v[i].z); }

    pub fn set(&mut self, i:usize, x: f32, y: f32, z: f32) {
        wrap_coordinate_to_box!(x, self.box_len, self.v[i].x);
        wrap_coordinate_to_box!(y, self.box_len, self.v[i].y);
        wrap_coordinate_to_box!(z, self.box_len, self.v[i].z);
    }

    pub fn add(&mut self, i:usize, x: f32, y: f32, z: f32) {
        wrap_coordinate_to_box!(self.v[i].x + x, self.box_len, self.v[i].x);
        wrap_coordinate_to_box!(self.v[i].y + y, self.box_len, self.v[i].y);
        wrap_coordinate_to_box!(self.v[i].z + z, self.box_len, self.v[i].z);
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
                                 i+1, " CA ", i+1, chain.x(i), chain.y(i), chain.z(i)).as_bytes()).ok();
    }
    out_writer.write(b"ENDMDL\n").ok();
}