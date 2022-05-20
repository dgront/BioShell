use rand::Rng;

use bioshell_ff::Coordinates;

pub fn random_chain(bond_length:f64, repulsion_distance: f64, system: &mut Coordinates) {

    let cutoff2 = (repulsion_distance * repulsion_distance) as f32;
    let (x, y, z) = random_unit_versor();
    system[1].x = system[0].x + (x * bond_length) as f32;
    system[1].y = system[0].y + (y * bond_length) as f32;
    system[1].z = system[0].z + (z * bond_length) as f32;

    for i in 2..system.size() {
        let mut go_on:bool = true;
        while go_on {
            let (x, y, z) = random_unit_versor();
            system[i].x = system[i-1].x + (x * bond_length) as f32;
            system[i].y = system[i-1].y + (y * bond_length) as f32;
            system[i].z = system[i-1].z + (z * bond_length) as f32;
            go_on = false;
            for j in 0..(i-1) {
                if system.distance_square(j,i) <= cutoff2 {
                    go_on = true;
                    break;
                }
            }
        }
    }
}


fn random_unit_versor() -> (f64, f64, f64) {

    let mut rng = rand::thread_rng();
    let x : f64 = rng.gen_range(-1.0..1.0);
    let y : f64 = rng.gen_range(-1.0..1.0);
    let z : f64 = rng.gen_range(-1.0..1.0);
    let l =  { (x * x + y * y + z * z).sqrt() };
    return (x/l, y/l, z/l);
}