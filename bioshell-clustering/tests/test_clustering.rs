extern crate core;

use std::iter::zip;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use nalgebra::{DMatrix, DVector};

use bioshell_clustering::{NeighborsOf, CartesianPoints};
use bioshell_clustering::optics::{Optics};
use bioshell_clustering::em::{expectation_maximization, log_likelihood};
use bioshell_clustering::kmeans::KMeans;
use bioshell_datastructures::euclidean_distance_squared;
use bioshell_statistics::{Distribution, MultiNormalDistribution, NormalDistribution};

#[test]
fn small_test_expectation_maximization() {

    let true_normals = vec![NormalDistribution::new(0.0, 1.0),
                           NormalDistribution::new(2.0, 0.5)];
    let true_weights = vec![1.0/3.0, 2.0/3.0];
    // --- data from true_normals
    let pts=[-0.538, 1.983, 1.305,-0.209, 2.139, 0.249,-0.896, 0.243,-1.240, 0.054, 1.750,
        2.402, 1.900, 1.247, 2.893, 2.239, 2.521, 2.379, 2.289, 1.696, 1.885, 1.775, 1.445, 2.605,
        2.519, 2.190, 0.987, 1.486, 2.233, 1.571];
    let mut data = vec![vec![0.0]; pts.len()];
    let true_loglikhd = log_likelihood(&true_normals, &true_weights, &data);

    for (i, p) in pts.iter().enumerate() { data[i][0] = *p; }

    // ---------- Distributions we will recover
    let mut normals = vec![NormalDistribution::new(0.5, 0.5),
                           NormalDistribution::new(1.8, 1.0)];

    // --- the true weight are: [1/3, 2/3]
    let mut weights = vec![0.5, 0.5];
    let loglikhd = expectation_maximization(&mut normals, &data, &mut weights, 1.0e-6);
    assert!(loglikhd >=true_loglikhd);
    // println!("{} {} {:?}", loglikhd, true_loglikhd, weights);
}

#[test]
fn large_test_expectation_maximization() {

    let n_dist: usize = 3;
    const N: usize = 1000;
    let n_data = [N*2, N*3, N*5];

    // ---------- Allocate space for training data
    let mut sample: Vec<Vec<f64>> = vec![vec!(0.0; 1); n_data[0] + n_data[1] + n_data[2]];
    // ---------- Distributions we will recover
    let mu = vec![0.0, 2.0, 5.0];
    let sigma = vec![0.5, 0.5, 0.5];
    let mut normals = vec![NormalDistribution::new(mu[0], sigma[0]),
                           NormalDistribution::new(mu[1], sigma[1]),
                           NormalDistribution::new(mu[2], sigma[2])];
    // ---------- Prepare training data
    let mut rng = SmallRng::seed_from_u64(0);
    let mut last: usize = 0;
    for i in 0..n_dist {
        for _j in 0..n_data[i] {
            normals[i].sample(&mut rng, &mut sample[last]);
            last += 1;
        }
    }
    let ref_likhd = log_likelihood(&normals, &vec![0.2, 0.3, 0.5], &sample);

    // --- change the distribution to make it harder
    for i in 0..n_dist {
        normals[i].set_parameters(mu[i] + 0.5, sigma[i]);
    }
    let mut weights = vec![0.0;n_dist];
    let log_lkhd = expectation_maximization(&mut normals, &sample, &mut weights, 1.0e-5);
    assert!(((log_lkhd-ref_likhd)/ref_likhd) < 0.01);
    // for (w,d) in zip(&weights, &normals) { println!("{} {}", w, d); }
}

#[allow(non_snake_case)]
#[test]
fn create_CartesianPoints() {

    const N: usize = 100;
    let mut points: Vec<[f64; 1]> = vec![[0.0]; N];
    for i in 0..N { points[i][0] = i as f64 * 0.1; }
    let ep = CartesianPoints::new(euclidean_distance_squared, points, 1);
    let nb = ep.neighbors_of(4,0.041);
    assert_eq!(nb.len(), 5);

    let mut nb_ids: Vec<usize> = nb.iter().map(|n| {n.idx}).collect();
    nb_ids.sort();
    let expected = vec![2, 3, 4, 5];
    for (actual, expcted) in zip(&nb_ids, &expected) {
        assert_eq!(actual, expcted);
    }
    for ni in nb_ids { println!("{}",ni)}
}

#[allow(non_snake_case)]
#[test]
fn Optics_clustering_Gaussian_data() {

    let mut n1: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu1 = DVector::<f64>::from_vec(vec![0.0, 0.0]);
    let sig1 = DMatrix::<f64>::from_iterator(2, 2, vec![0.1, 0.05, 0.05, 0.1].into_iter());
    n1.set_parameters(&mu1, &sig1);

    let mut n2: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu2 = DVector::<f64>::from_vec(vec![3.0, 3.0]);
    let sig2 = DMatrix::<f64>::from_iterator(2, 2, vec![0.1, 0.03, 0.03, 0.1].into_iter());
    n2.set_parameters(&mu2, &sig2);

    // --- container for the data to be clustered
    let mut data: Vec<Vec<f64>> = Vec::new();

    let mut rng = SmallRng::seed_from_u64(0);
    for _ in 0..20 {
        let mut row = vec![0.0, 0.0];
        n1.sample(&mut rng, &mut row);
        data.push(row);
        row = vec![0.0, 0.0];
        n2.sample(&mut rng, &mut row);
        data.push(row);
    }

    let opt_clust = Optics::new(0.5, 5,
                    Box::new(CartesianPoints::new(euclidean_distance_squared, data.clone(), 2)));

    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    assert_eq!(cluster_size[0], 20);
    assert_eq!(cluster_size[1], 20);
}

#[allow(non_snake_case)]
#[test]
fn test_Optics() {
    let data: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0],
                                   vec![2.0, 0.0], vec![3.0, 0.0], vec![4.0, 0.0], vec![6.0, 0.0],
                                   vec![8.0, 0.0], vec![10.0, 0.0], vec![10.0, 1.0], vec![10.0, -1.0],
                                   vec![18.0, 0.0]];


    let opt_clust = Optics::new(2.5, 4,
        Box::new(CartesianPoints::new(euclidean_distance_squared, data.clone(), 2)));
    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    assert_eq!(data.len(), cluster_size.iter().sum());
    let expected_size: Vec<usize> = vec![8, 4, 1];
    let expected_order: Vec<usize> = vec![0, 1, 3, 4, 2, 5, 6, 7, 8, 9, 11, 10, 12];
    assert_eq!(expected_size, *cluster_size);
    assert_eq!(expected_order, *opt_clust.clustering_order());
}

#[test]
fn test_kmeans() {
    let data: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![0.0, 1.0], vec![1.0, 0.0], vec![1.0, 1.0], vec![1.0, 2.0],
                                   vec![2.0, 1.0], vec![3.0, 3.0], vec![4.0, 2.0], vec![4.0, 3.0],
                                   vec![4.0, 3.5], vec![4.0, 5.0], vec![5.0, 4.0], vec![5.0, 5.0],];
    let mut kmeans = KMeans::new(2, data, 2, euclidean_distance_squared);
    kmeans.cluster_n(0.01, 100);
    let mut out = kmeans.assignments().clone();
    if out[0] == 0 { out.reverse();}
    for i in 0..6 { assert_eq!(out[i], 1); }
    for i in 7..out.len() { assert_eq!(out[i], 0); }
}