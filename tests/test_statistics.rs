use bioshell_numerical::statistics::{Distribution, Histogram, MultiNormalDistribution, OnlineMultivariateStatistics};
use rand_distr::{Normal, Distribution as RndDistribution};
use nalgebra::{DMatrix, DVector};

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

#[test]
fn test_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    let mu = DVector::<f64>::repeat(2, 0.1);
    let sig = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.1, 0.1, 2.0].into_iter());
    n.set_parameters(&mu, &sig);
    let logprob = n.logpdf(&vec![1.0, 0.0]);
    assert!((logprob + 3.469636899044226).abs() < 0.0001);
}

#[test]
fn sample_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
                     &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));
    let n_samples = 100000;
    let n_dim: usize = 2;

    let mut stats = OnlineMultivariateStatistics::new(n_dim);
    let mut row = vec!(0.0; n_dim);

    for _ in 0..n_samples {
        n.rand(&mut row);
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


#[test]
fn estimate_MultiNormalDistribution() {
    let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
                     &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));

    let mu = DVector::<f64>::repeat(2, 0.1);
    let sig = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.1, 0.1, 2.0].into_iter());
    n.set_parameters(&mu, &sig);
    let logprob = n.logpdf(&vec![1.0, 0.0]);
    assert!((logprob + 3.469636899044226).abs() < 0.0001);
}

#[test]
fn test_OnlineMultivariateStatistics() {
    let n_samples = 100000;
    let n_dim: usize = 4;

    let mut normal = Normal::new(2.0, 3.0).unwrap();         // --- mean 2, standard deviation 3
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