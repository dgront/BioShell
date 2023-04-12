use std::time::Instant;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use bioshell_clustering::kd_tree::{count, create_kd_tree, euclidean_distance_squared};

fn construct_k3_tree() {

    const N: usize = (2_usize.pow(20) - 1) as usize;
    let mut rng = SmallRng::seed_from_u64(0);
    let mut data = vec![vec![0.0]; N];
    for i in 0..data.len() { data[i][0] = rng.gen(); }

    let query = vec![0.5];
    // --- find nearest with brute force
    let (mut min_d, mut min_e) = (euclidean_distance_squared(&query, &data[0], 1), &data[0]);
    for e in data.iter() {
        let d = euclidean_distance_squared(&query, e, 1);
        if d < min_d { (min_d, min_e) = (d, e);}
    }

    let start = Instant::now();
    let root = create_kd_tree(&mut data.clone(), 1).unwrap();
    let end = start.elapsed();

    println!("construct_k3_tree(): {:.2?}", end);
}

fn main() {
    construct_k3_tree();
}