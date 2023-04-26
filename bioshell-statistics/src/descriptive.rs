/// Provides on-line statistics for N-dimensional samples
///
/// This struct accumulates observations without actually storing them, and on the fly provides
/// basic descriptive parameters of the accumulated sample
///
/// # Examples
/// The following example:
/// - creates a 2-dimensional normal distribution
/// - creates an `OnlineMultivariateStatistics` object
/// - withdraws 50000 observations from the distribution and accumulates them
/// - checks the convergence of expectations vector
///
/// Refer to [`MultiNormalDistribution`](MultiNormalDistribution) struct documentation to see
/// how to initialize parameters of the distribution.
/// ```
/// use bioshell_statistics::{OnlineMultivariateStatistics, MultiNormalDistribution, Distribution};
/// # use nalgebra::{DMatrix, DVector};
/// # use rand::rngs::SmallRng;
/// # use rand::SeedableRng;
///
/// # let mu = DVector::from_vec(vec![1.0, 2.0]);
/// # let sig = DMatrix::<f64>::from_iterator(2, 2, vec![0.2, 0.1, 0.1, 2.0].into_iter());
/// let mut gen: MultiNormalDistribution = MultiNormalDistribution::new(2);
/// // ------ initialize MultiNormalDistribution here ...
/// # gen.set_parameters(&mu, &sig);
/// # let mut stats = OnlineMultivariateStatistics::new(2);
/// let mut observation = [0.0, 0.0];
/// let mut rng = SmallRng::seed_from_u64(0);
/// for i in 0..50000 {
///         gen.sample(&mut rng, &mut observation);   // --- fill observation vector with new random values
///         stats.accumulate(&observation);             // --- accumulate them
///         if i % 1000 == 0 {
///             println!("Average converged to {:.4}, {:.4}", stats.avg()[0], stats.avg()[1]);
///         }
/// }
/// assert!((stats.avg()[0] - 1.0).abs() < 0.01);
/// assert!((stats.avg()[1] - 2.0).abs() < 0.01);
/// ```
#[derive(Clone)]
pub struct OnlineMultivariateStatistics {
    dim: usize,
    count: usize,
    m1: Vec<f64>,
    m2: Vec<f64>,
    min: Vec<f64>,
    max: Vec<f64>,
    cov: Vec<Vec<f64>>,
}

impl OnlineMultivariateStatistics {

    /// Create a new object to gather statistics on N-dimensional samples
    ///
    /// # Arguments
    /// * `dim` - dimensionality of the data being observed
    ///
    pub fn new(dim: usize) -> OnlineMultivariateStatistics {
        OnlineMultivariateStatistics{dim, count:0, m1: vec![0.0; dim], m2: vec![0.0; dim],
            min: vec![0.0; dim], max: vec![0.0; dim],
            cov: vec![vec![0.0; dim]; dim]}
    }

    /// Returns the number of dimensions of the observed vectors
    pub fn dim(&self) -> usize { self.dim }

    /// Accumulate a single N-dimensional point
    pub fn accumulate(&mut self, d:&[f64]) {

        assert_eq!(d.len(), self.dim);                  // --- incoming vector must be of the same size at the statistics

        if self.count == 0 {                            // --- if this is the very first observation...
            for i in 0..self.dim {
                self.min[i] = d[i];                     // --- copy it as min and max
                self.max[i] = d[i];
            }
        }
        self.count += 1;
        for i in 0..self.dim {
            let delta_x: f64 = d[i] - self.m1[i];
            self.m1[i] += delta_x / self.count as f64;  // --- M1[i] is now the new average for i-th dimension
            self.m2[i] += delta_x * (d[i] -self.m1[i]);

            for j in i+1..self.dim {             // --- here we update the covariance terms
                let delta_y: f64 = d[j] - self.m1[j];
                self.cov[i][j] += delta_y * (d[i] - self.m1[i]);
                self.cov[j][i] = self.cov[i][j];
            }

            self.min[i] = self.min[i].min(d[i]);
            self.max[i] = self.max[i].max(d[i]);
        }
    }

    /// Accumulate a single 1-dimensional point
    pub fn accumulate_1d(&mut self, x:f64) {
        let v = [x];
        self.accumulate(&v);
    }

    /// Returns the number of observed samples
    pub fn count(&self) -> usize{ self.count }

    /// Returns the vector of minimum values observed in each of dimensions
    ///
    /// Note, that the maximum values were taken independently for each dimension and the reported values
    /// come from separate observations.
    pub fn min(&self) -> &Vec<f64> { &self.min }

    /// Returns the vector of maximum values observed in each of dimensions
    ///
    /// Note, that the maximum values were taken independently for each dimension and the reported values
    /// come from separate observations.
    pub fn max(&self) -> &Vec<f64> { &self.max }

    /// Returns the average vector of observations
    pub fn avg(&self) -> &Vec<f64> { &self.m1 }

    /// Returns the vector of variance values for each dimension
    pub fn var(&self) -> Vec<f64> {
        let mut ret = self.m2.clone();
        for i in 0..ret.len() {
            ret[i] /= self.count as f64 - 1.0;
        }

        return ret;
    }

    /// Returns the covariance matrix
    pub fn cov(&self) -> Vec<Vec<f64>> {

        let mut ret = self.cov.clone();
        for i in 0..ret.len() {
            for j in 0..i {
                ret[i][j] /= self.count as f64 - 1.0;
                ret[j][i] /= self.count as f64 - 1.0;
            }
            ret[i][i] = self.m2[i] / (self.count as f64 - 1.0);
        }

        return ret;
    }
}