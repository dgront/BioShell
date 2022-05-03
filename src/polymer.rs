use rand::Rng;
use std::env;
use std::time::Instant;

use bioshell_numerical::Vec3;

pub fn main() {
    let args: Vec<String> = env::args().collect();

    let v = Vec3::new(0.0,0.0,0.0);
    println!("{:.2} {:.2} {:.2}", v.x, v.y, v.z );
}
