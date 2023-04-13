use std::time::Instant;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use bioshell_clustering::kd_tree::{create_kd_tree, find_within};
use bioshell_clustering::{euclidean_distance_squared};

fn construct_k3_tree() {

    const N_PTS: usize = (2_usize.pow(20) - 1) as usize;
    const N_DIM: usize = 3;
    const N_SEARCH: usize = 10000;

    let mut rng = SmallRng::seed_from_u64(0);
    let mut data = vec![vec![0.0; N_DIM]; N_PTS];
    for i in 0..data.len() {
        for j in 0..N_DIM { data[i][j] = rng.gen(); }
    }

    let start = Instant::now();
    let root = create_kd_tree(&mut data.clone(), 1).unwrap();
    println!("construct_k3_tree() of {} points: {:.2?}", N_PTS, start.elapsed());

    // --- check how quickly can we find neighbors; note that the search radius is the squared value of the euclidean distance
    let start = Instant::now();
    let radius = 0.01;
    for i in 0..N_SEARCH {
        find_within(&root, &data[i], N_DIM, radius * radius, euclidean_distance_squared);
    }
    println!("find_within() x {}: {:.2?}", N_SEARCH, start.elapsed());
}

fn main() {
    construct_k3_tree();
}