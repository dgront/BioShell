//! K-means clustering algorithm splits observations into separate K groups.
//!
use std::ops::{IndexMut};
use log::debug;
use rand::{Rng, SeedableRng};
use rand::rngs::SmallRng;

/// Provides the k-means clustering algorithm.
///
/// K-means method partitions a set of unlabeled observations into clusters,
/// where each observation belongs to the cluster with the nearest mean. As a result, the code
/// provides assignment of each point to one of k-clusters. The value of `k` must be provided by a user.
/// In the current implementation provides only random initialisation, i.e. the `k` starting points
/// are randomly selected from the data set. Since the final clustering result is highly affected by that
/// initial assignment, it's advised to repeat the procedure to find the best solution.
pub struct KMeans<T, D> {
    /// The total number of points to be clustered
    n: usize,
    /// The number of clusters to be created
    k: usize,
    /// number of dimensions for each point of the type T
    dimensionality: usize,
    /// Data structure that holds the input points - vector of size N, that holds N input points
    points: Vec<T>,
    /// Cluster centers - vector of size K
    centers: Vec<T>,
    /// holds the shortest distance between each point and it's cluster center - vector of size N
    closest_dist: Vec<f64>,
    /// holds the index of the cluster a given i-th points belongs to - vector of size N
    closest_center: Vec<usize>,
    /// holds the size of each cluster - vector of size K
    cluster_counts: Vec<usize>,
    /// the distance function
    distance: D
}

impl<T, D> KMeans<T, D> where T: IndexMut<usize, Output = f64> + Clone, D: Fn(&T, &T, usize) -> f64 {

    pub fn new(k:usize, data:Vec<T>, dimensionality: usize, distance: D) -> KMeans<T, D> {

        let centers = vec![data[0].clone(); k];
        let n = data.len();
        KMeans{
            n, k, dimensionality,
            points: data,
            centers,
            closest_dist: vec![0.0; n],
            closest_center: vec![0; n],
            cluster_counts: vec![0; k],
            distance
        }
    }

    /// Run the clustering procedure
    ///
    /// # Arguments
    /// * `epsilon` - convergence criterion: a clustering stops when the relative change in error function
    ///     is smaller that `epsilon`. The error function is defined as the sum of distances from each
    ///     point to its respective cluster center.
    ///
    /// # Returns
    ///     error value reached in the last iteration
    pub fn cluster(&mut self, epsilon: f64) -> f64 {
        // --- initialize by assigning k initial clusters
        // self.init_random();
        // --- initialize according to the kmeans++ method
        self.init_plus_plus();
        // --- initial assignment points to clusters
        let mut new_err = self.assign();
        let mut prev_err: f64 = new_err * 2.0;
        while (prev_err - new_err) / prev_err > epsilon {
            prev_err = new_err;
            // --- compute the mean for each cluster
            for c in 0..self.k {
                self.cluster_counts[c] = 0;
                for d in 0..self.dimensionality { self.centers[c][d] = 0.0; }
            }
            for i in 0..self.n {
                let which_center = self.closest_center[i];
                self.cluster_counts[which_center] += 1;
                for d in 0..self.dimensionality { self.centers[which_center][d] += self.points[i][d]; }
            }
            for c in 0..self.k {
                for d in 0..self.dimensionality { self.centers[c][d] /= self.cluster_counts[c] as f64; }
            }
            // --- re-assign points
            new_err = self.assign();
        }
        return new_err;
    }

    /// Repeat the [`cluster()`](cluster()) call multiple times to find a better solution
    ///
    /// Repeats the clustering `n_repeats` times and records the best cluster assignment found. The lowest
    /// value of the sum of distances is returned
    pub fn cluster_n(&mut self, epsilon: f64, n_repeats: usize) -> f64 {

        let mut err = self.cluster(epsilon);
        let mut assignmnt = self.assignments().clone();
        let mut centers = self.centers.clone();
        for i in 1..n_repeats {
            let new_err = self.cluster(epsilon);
            debug!("min-distance after {} iteration: {}", i, err);
            if new_err < err {
                err = new_err;
                assignmnt = self.assignments().clone();
                centers = self.centers.clone();
            }
        }
        self.closest_center = assignmnt;
        self.centers = centers;

        return err;
    }

    /// Provides assignments of points to clusters
    ///
    /// Returned vector of size `N` provides
    pub fn assignments(&self) -> &Vec<usize> { &self.closest_center }

    /// Provides centers of the `k` clusters found in the most recent clustering run.
    ///
    /// Each cluster center is represented by a variable of type `T` - the same as the input points
    pub fn centers(&self) -> &Vec<T> { &self.centers }

    pub fn sizes(&self) -> Vec<usize> {
        let mut ret: Vec<usize> = vec![0; self.k];
        for a in &self.closest_center { ret[*a] += 1; }
        return ret;
    }

    fn init_random(&mut self) {
        let mut rng = SmallRng::from_entropy();
        // let mut rng = SmallRng::seed_from_u64(115);

        let mut selected_pts = vec![];
        let mut which = rng.gen_range(0..self.n);
        selected_pts.push(which);
        for i in 1..self.k {
            which = rng.gen_range(0..self.n);
            while selected_pts.contains(&which) {
                which = rng.gen_range(0..self.n);
            }
            selected_pts.push(which);
            self.centers[i] = self.points[which].clone();
        }
    }

    fn init_plus_plus(&mut self) {
        let mut rng = SmallRng::from_entropy();

        // ------ select the first cluster center randomly and push it to the list
        let which = rng.gen_range(0..self.n);
        self.centers[0] = self.points[which].clone();

        let mut distances: Vec<f64> = vec![0.0; self.n];
        for k in 1..self.k {                // --- find k-1 additional cluster centers
            let mut sum_d = 0.0;
            for i in 0..self.n {            // --- minimum distance between i-th point and any center
                let mut min_d = (self.distance)(&self.centers[0], &self.points[i], self.dimensionality);
                for j in 1..k {
                    min_d = (self.distance)(&self.centers[j], &self.points[i], self.dimensionality).min(min_d);
                }
                distances[i] = min_d * min_d;
                sum_d += min_d * min_d;
            }
            // ------ select the next cluster by the kmeans++ rule, i.e. proportional to D(x)^2
            let w = rng.gen_range(0.0..sum_d);
            let mut k_which = 0;
            let mut k_sum = distances[0];
            while k_sum < w {
                k_which += 1;
                k_sum += distances[k_which];
            }
            if k_which==self.k { k_which -= 1; }
            self.centers[k] = self.points[k_which].clone();
        }
    }

    fn assign(&mut self) -> f64 {
        for i in 0..self.n {         // --- initialise distances by the first cluster
            self.closest_dist[i] = (self.distance)(&self.centers[0], &self.points[i], self.dimensionality);
            self.closest_center[i] = 0;
        }
        for j in 1..self.k {             // --- iterate over clusters
            for i in 0..self.n {         // --- iterate over points
                let d: f64 = (self.distance)(&self.centers[j], &self.points[i], self.dimensionality);
                if self.closest_dist[i] > d {
                    self.closest_dist[i] = d;
                    self.closest_center[i] = j;
                }
            }
        }
        self.closest_dist.iter().sum()
    }
}