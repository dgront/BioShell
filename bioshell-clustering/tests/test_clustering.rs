use rand_distr::{Normal, Distribution as RndDistribution};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use nalgebra::{DMatrix, DVector};
use rand::Rng;

use bioshell_clustering::{expectation_maximization, Optics, EuclideanPoints};

use bioshell_statistics::{Distribution, MultiNormalDistribution,
                                      NormalDistribution, OnlineMultivariateStatistics};


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
                    Box::new(EuclideanPoints::new(data.clone())));

    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    println!("{:?}",cluster_size);
}

#[allow(non_snake_case)]
#[test]
fn test_Optics() {
    let data: Vec<Vec<f64>> = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0],
                                   vec![2.0, 0.0], vec![3.0, 0.0], vec![4.0, 0.0], vec![6.0, 0.0],
                                   vec![8.0, 0.0], vec![10.0, 0.0], vec![10.0, 1.0], vec![10.0, -1.0],
                                   vec![18.0, 0.0]];


    let opt_clust = Optics::new(2.5, 4,
                                    Box::new(EuclideanPoints::new(data.clone())));
    let cluster_size: Vec<usize> = opt_clust.clusters().iter().map(|r| r.len()).collect();
    assert_eq!(data.len(), cluster_size.iter().sum());
    let expected_size: Vec<usize> = vec![8, 4, 1];
    let expected_order: Vec<usize> = vec![0, 1, 3, 4, 2, 5, 6, 7, 8, 9, 11, 10, 12];
    assert_eq!(expected_size, *cluster_size);
    assert_eq!(expected_order, *opt_clust.clustering_order());
}