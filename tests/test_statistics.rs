use rand_distr::{Normal, Distribution as RndDistribution};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use nalgebra::{DMatrix, DVector};
use rand::Rng;
use std::string::String;
use std::fmt::Write;

use bioshell_numerical::statistics::{Distribution, Estimable, Histogram, MultiNormalDistribution,
        NormalDistribution, OnlineMultivariateStatistics, expectation_maximization};

use bioshell_numerical::clustering::{Euclidean, Optics};

/// Create a histogram and fill it with deterministic observations
#[test]
fn create_histogram() {
    let test_data = vec!(1.0, 1.1, 1.3, 1.6, 1.7, 2.0);
    let mut h: Histogram = Histogram::by_bin_width(0.5);
    for x in test_data { h.insert(x); }
    assert_eq!(h.which_bin(1.11), 2);
    assert_eq!(h.which_bin(1.49), 2);
    assert_eq!(h.which_bin(1.51), 3);
    assert_eq!(h.sum(), 6.0);
}

/// Check if MultiNormalDistribution correctly evaluates its probability
#[test]
#[allow(non_snake_case)]
fn test_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu = DVector::<f64>::repeat(2, 0.1);
    let sig = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.1, 0.1, 2.0].into_iter());
    n.set_parameters(&mu, &sig);
    let logprob = n.logpdf(&vec![1.0, 0.0]);
    assert!((logprob + 3.469636899044226).abs() < 0.0001);
}

/// sample from MultiNormalDistribution, check is it correctly estimates its parameters
#[allow(non_snake_case)]
#[test]
fn sample_MultiNormalDistribution() {

    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
                     &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));
    let n_samples = 100000;
    let n_dim: usize = 2;

    let mut stats = OnlineMultivariateStatistics::new(n_dim);
    let mut row = vec!(0.0; n_dim);

    let mut rng = SmallRng::seed_from_u64(0);
    for _ in 0..n_samples {
        n.sample(&mut rng, &mut row);
        stats.accumulate(&row);
    }
    println!("{} {}",stats.covar(0,0) ,stats.covar(0,1) );
    assert!((stats.avg(0) - 1.0).abs() < 0.01);
    assert!((stats.avg(1) - 2.0).abs() < 0.01);
    assert!((stats.covar(0,1) - 0.5).abs() < 0.01);
    assert!((stats.covar(1,0) - 0.5).abs() < 0.01);
    assert!((stats.var(0).sqrt() - 1.0).abs() < 0.01);
    assert!((stats.var(1).sqrt() - 1.0).abs() < 0.01);
}

/// Print a MultiNormalDistribution object using fmt::Display
#[allow(non_snake_case)]
#[test]
fn format_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
                     &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));

    let expected = "mu =  [ 1.0000,  2.0000], sigma = [ [ 1.0000,  0.5000], [ 0.5000,  1.0000]]";
    let mut actual = String::new();
    write!(actual, "{}", n).unwrap();
    assert_eq!(actual, expected);
}

#[allow(non_snake_case)]
#[test]
fn test_NormalDistribution() {
    // --- create a N(2.0, 1.5) distribution
    let mut n: NormalDistribution = NormalDistribution::new(2.0, 1.5);

    // --- evaluate pdf
    let mut prob = n.pdf(&vec![1.0]);
    assert!((prob - 0.21297).abs() < 0.0001);

    // --- reset the distribution and evaluate pdf again
    n.set_parameters(0.0, 3.0);
    prob = n.pdf(&vec![1.0]);
    assert!((prob - 0.12579).abs() < 0.0001);
    assert!((n.mean()).abs() < 0.0001);
    assert!((n.sdev() - 3.0).abs() < 0.0001);

    // --- estimate the distribution from a sample
    let n_samples = 1000000;
    let mut sample: Vec<Vec<f64>> = vec!(vec!(0.0; 1); n_samples);
    let mut rng = SmallRng::seed_from_u64(0);
    for v in sample.iter_mut() { n.sample(&mut rng,v); }
    n.estimate(&sample);
    assert!((n.mean()).abs() < 0.1);
    assert!((n.sdev() - 3.0).abs() < 0.1);
}

#[test]
fn test_expectation_maximization() {

    let n_dist: usize = 3;
    let n_data_1: usize = 1000;
    let n_data: usize = n_data_1 * n_dist;

    // ---------- Allocate space for training data
    let mut sample: Vec<Vec<f64>> = vec!(vec!(0.0; 1); n_data);
    // ---------- Distributions we will recover
    let mut normals = vec![NormalDistribution::new(0.0, 0.5),
                           NormalDistribution::new(2.0, 0.5),
                           NormalDistribution::new(5.0, 0.5)];
    // ---------- Prepare training data
    let mut rng = SmallRng::seed_from_u64(0);
    for i in 0..n_dist {
        for j in 0..n_data_1 { normals[i].sample(&mut rng, &mut sample[i * n_data_1 + j]); }
        normals[i].set_parameters(2.5, 3.0);    // --- change the distribution to make it harder
    }
    // ---------- Randomly assign data points to distributions
    let mut assignment: Vec<Vec<bool>> = vec!(vec!(false; n_data); n_dist);
    for i in 0..n_data {
        let i_dist = rand::thread_rng().gen_range(0..n_dist);
        assignment[i_dist][i] = true;
    }

    expectation_maximization(&mut normals, &sample, &mut assignment, 0.000001);
    let mut totals: Vec<i32> = vec!(0; n_dist);
    for i in 0..n_dist {
        totals[i] = (&assignment[i]).into_iter().map(|v| *v as i32).sum();
    }
    println!("{:?}",&totals);
}

#[allow(non_snake_case)]
#[test]
fn test_OnlineMultivariateStatistics() {
    let n_samples = 100000;
    let n_dim: usize = 4;

    let normal = Normal::new(2.0, 3.0).unwrap();         // --- mean 2, standard deviation 3
    let mut stats = OnlineMultivariateStatistics::new(n_dim);
    let mut row = vec!(0.0; n_dim);
    for _ in 0..n_samples {
        for i in 0..row.len() {
            row[i] = normal.sample(&mut rand::thread_rng());
        }
        stats.accumulate(&row);
    }
    for i in 0..n_dim {
        assert!((stats.avg(i)-2.0).abs() < 0.1);
        assert!((stats.var(i).sqrt()-3.0).abs() < 0.1);
    }
}

#[allow(non_snake_case)]
#[test]
fn Optics_clustering_Gaussian_data() {

    let mut n1: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu1 = DVector::<f64>::from_vec(vec![0.0, 0.0]);
    let sig1 = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.1, 0.1, 0.2].into_iter());
    n1.set_parameters(&mu1, &sig1);

    let mut n2: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu2 = DVector::<f64>::from_vec(vec![3.0, 3.0]);
    let sig2 = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.03, 0.03, 0.2].into_iter());
    n2.set_parameters(&mu2, &sig2);

    // --- container for the data to be clustered
    let mut data: Vec<Vec<f64>> = Vec::new();

    let mut rng = SmallRng::seed_from_u64(0);
    for _ in 0..100 {
        let mut row = vec![0.0, 0.0];
        n1.sample(&mut rng, &mut row);
        data.push(row);
        row = vec![0.0, 0.0];
        n2.sample(&mut rng, &mut row);
        data.push(row);
    }

    let mut opt_clust = Optics::new(0.5, 5, Euclidean{});
    opt_clust.run_clustering(&data);

    let &d =  &opt_clust.reacheability_distance();
    let &o =  &opt_clust.clustering_order();
    let mut cluster_index = 0;
    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    for ci in opt_clust.clusters() {
        for i in ci.clone() {
            let ei:usize = o[i];
            println!("c: {} {} {} {}  {}", data[ei][0], data[ei][1], cluster_index, d[i], ei);
        }
        cluster_index += 1;
    }
    println!("{:?}",cluster_size);
}

#[allow(non_snake_case)]
#[test]
fn test_Optics() {
    let data: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0],
                                   vec![2.0, 0.0], vec![3.0, 0.0], vec![4.0, 0.0], vec![6.0, 0.0],
                                   vec![8.0, 0.0], vec![10.0, 0.0], vec![10.0, 1.0], vec![10.0, -1.0],
                                   vec![18.0, 0.0]];


    let mut opt_clust = Optics::new(2.5, 4, Euclidean{});
    opt_clust.run_clustering(&data);
    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    assert_eq!(data.len(), cluster_size.iter().sum());
    let expected_size: Vec<usize> = vec![8, 4, 1];
    let expected_order: Vec<usize> = vec![0, 1, 3, 4, 2, 5, 6, 7, 8, 9, 11, 10, 12];
    assert_eq!(expected_size, *cluster_size);
    assert_eq!(expected_order, *opt_clust.clustering_order());
}