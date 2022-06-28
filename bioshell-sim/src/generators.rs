use rand::Rng;

use bioshell_ff::{Coordinates};
use bioshell_numerical::Vec3;

pub fn random_chain(bond_length:f64, repulsion_distance: f64, start: &Vec3, system: &mut Coordinates) {

    let cutoff2 = (repulsion_distance * repulsion_distance) as f32;

    system.set_x(0, start.x);
    system.set_y(0, start.y);
    system.set_z(0, start.z);

    let (x, y, z) = random_unit_versor();
    system.set_x(1, system.x(0) + (x * bond_length) as f32);
    system.set_y(1, system.y(0) + (y * bond_length) as f32);
    system.set_z(1, system.z(0) + (z * bond_length) as f32);

    for i in 2..system.size() {
        let mut go_on: bool = true;
        while go_on {
            let (x, y, z) = random_unit_versor();
            system.set_x(i, system.x(i-1) + (x * bond_length) as f32);
            system.set_y(i, system.y(i-1) + (y * bond_length) as f32);
            system.set_z(i, system.z(i-1) + (z * bond_length) as f32);
            go_on = false;
            for j in 0..(i-1) {
                if system.closest_distance_square(j, i) <= cutoff2 {
                    go_on = true;
                    break;
                }
            }
        }
    }
}

pub fn cubic_grid_atoms(system: &mut Coordinates) {

    let points_one_side: usize = (f64::powf(system.size() as f64, 1.0 / 3.0)).ceil() as usize;
    let dw = system.box_len() / points_one_side as f32;
    let cell_margin = dw / 2.0;

    for i in 0..system.size() {
        let k = i % points_one_side;
        let l = (i / points_one_side) % points_one_side;
        let m = (i / (points_one_side * points_one_side)) % points_one_side;
        system.set(i,dw * k as f32 + cell_margin,dw * l as f32 + cell_margin,dw * m as f32 + cell_margin)
    }
}


pub fn square_grid_atoms(system: &mut Coordinates) {

    let points_one_side: usize = (f64::powf(system.size() as f64, 0.5)).ceil() as usize;
    let dw = system.box_len() / points_one_side as f32;
    let cell_margin = dw / 2.0;

    for i in 0..system.size() {
        let k = i % points_one_side;
        let l = i / points_one_side;
        system.set(i,dw * k as f32 + cell_margin,dw * l as f32 + cell_margin,0.0);
    }
}

fn random_unit_versor() -> (f64, f64, f64) {

    let mut rng = rand::thread_rng();
    let x : f64 = rng.gen_range(-1.0..1.0);
    let y : f64 = rng.gen_range(-1.0..1.0);
    let z : f64 = rng.gen_range(-1.0..1.0);
    let l =  { (x * x + y * y + z * z).sqrt() };
    return ((x/l) as f64, (y/l) as f64, (z/l) as f64);
}

fn random_unit_versor32() -> (f32, f32, f32) {

    let mut rng = rand::thread_rng();
    let x : f64 = rng.gen_range(-1.0..1.0);
    let y : f64 = rng.gen_range(-1.0..1.0);
    let z : f64 = rng.gen_range(-1.0..1.0);
    let l =  { (x * x + y * y + z * z).sqrt() };
    return ((x/l) as f32, (y/l) as f32, (z/l) as f32);
}