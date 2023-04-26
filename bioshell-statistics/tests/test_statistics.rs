use rand_distr::{Normal, Distribution as RndDistribution};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use nalgebra::{DMatrix, DVector};
use std::string::String;
use std::fmt::Write;

use bioshell_statistics::{Distribution, Estimable, Histogram, MultiNormalDistribution,
        NormalDistribution, OnlineMultivariateStatistics};

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
    let mu = DVector::<f64>::from_vec(vec![0.1, 0.1]);
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

    let mut stats = OnlineMultivariateStatistics::new(2);
    let mut row = [0.0, 0.0];

    let mut rng = SmallRng::seed_from_u64(0);
    for _ in 0..n_samples {
        n.sample(&mut rng, &mut row);
        stats.accumulate(&row);
    }
    assert!((stats.avg()[0] - 1.0).abs() < 0.01);
    assert!((stats.avg()[1] - 2.0).abs() < 0.01);
    assert!((stats.cov()[0][1] - 0.5).abs() < 0.01);
    assert!((stats.cov()[1][0] - 0.5).abs() < 0.01);
    assert!((stats.var()[0].sqrt() - 1.0).abs() < 0.01);
    assert!((stats.var()[1].sqrt() - 1.0).abs() < 0.01);
}

/// Print a MultiNormalDistribution object using fmt::Display
#[allow(non_snake_case)]
#[test]
fn format_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
                     &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));

    let expected = "'mu': [ 1.0000,  2.0000], 'sigma': [ [ 1.0000,  0.5000], [ 0.5000,  1.0000]]";
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
        assert!((stats.avg()[i]-2.0).abs() < 0.1);
        assert!((stats.var()[i].sqrt()-3.0).abs() < 0.1);
    }
}

