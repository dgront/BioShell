#[derive(Clone)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub res_type: i8,
    pub atom_type: i8,
    pub chain_id: i16
}

impl Vec3 {
    pub fn new(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3 { x: x, y: y, z: z, res_type:0, atom_type: 0, chain_id: 0}
    }
    pub fn from_float(value: f32) -> Vec3 {
        Vec3 {
            x: value,
            y: value,
            z: value,
            res_type:0, atom_type: 0, chain_id: 0
        }
    }
}

fn main() {
    let v1 = Vec3::from_float(0.0);
    println!("{:.2} {:.2}", v1.x, v1.y );
}