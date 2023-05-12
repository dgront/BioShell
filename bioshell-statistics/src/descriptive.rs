use std::iter::zip;

/// Calculates the average value of a given sample
///
/// # Examples
/// ```
/// use bioshell_statistics::avg;
/// let a = avg(&[8.0, 9.0, 8.5]);
/// assert!((a - 8.5).abs() < 0.00001);
/// ```
pub fn avg(sample: &[f64]) -> f64 {
    let sum: f64 = sample.iter().sum();
    return sum / sample.len() as f64;
}

/// Calculates the weighted average value of a given sample
///
/// # Examples
/// ```
/// use bioshell_statistics::avg_weighted;
/// let a = avg_weighted(&[8.0, 9.0, 8.5], &[0.2, 0.3, 0.5]);
/// assert!((a - 8.55).abs() < 0.00001);
/// ```
pub fn avg_weighted(sample: &[f64], weights: &[f64]) -> f64 {
    let mut sum = 0.0;
    for (x, w) in zip(sample, weights) { sum += x * w; }
    return sum / weights.iter().sum::<f64>();
}

/// Calculates the weighted variance value of a given sample
///
/// # Examples
/// ```
/// use bioshell_statistics::var;
/// let v = var(&[3.0, 6.0, 5.0, 4.0, 7.0]);
/// assert!((v - 2.0).abs() < 0.00001);
/// ```
pub fn var(sample: &[f64]) -> f64 {

    let avg: f64 = avg(sample);
    let mut ret: f64 = 0.0;
    for x in sample { ret += (x - avg) * (x - avg); }
    return ret / sample.len() as f64;
}

/// Calculates the weighted variance value of a given sample
///
/// # Examples
/// ```
/// use bioshell_statistics::var_weighted;
/// let v = var_weighted(&[3.0, 6.0, 5.0, 4.0, 7.0], &[0.2, 0.3, 0.4, 0.5, 0.6]);
/// assert!((v - 1.91).abs() < 0.00001);
/// ```
pub fn var_weighted(sample: &[f64], weights: &[f64]) -> f64 {

    let avg: f64 = avg_weighted(sample, weights);
    let mut var: f64 = 0.0;
    for (x, w) in zip(sample, weights) { var += w * (x - avg) * (x - avg); }
    return var / weights.iter().sum::<f64>();
}

/// Calculates the weighted average value and the variance of a given sample
///
/// # Examples
/// ```
/// use bioshell_statistics::avg_var_weighted;
/// let (a, v) = avg_var_weighted(&[3.0, 6.0, 5.0, 4.0, 7.0], &[0.2, 0.3, 0.4, 0.5, 0.6]);
/// assert!((a - 5.3).abs() < 0.00001);
/// assert!((v - 1.91).abs() < 0.00001);
/// ```
pub fn avg_var_weighted(sample: &[f64], weights: &[f64]) -> (f64,f64) {

    let avg: f64 = avg_weighted(sample, weights);
    let mut var: f64 = 0.0;
    for (x, w) in zip(sample, weights) { var += w * (x - avg) * (x - avg); }
    return (avg,var / weights.iter().sum::<f64>());
}

/// Calculate covariance from two columns of data
pub fn cov(x: &[f64], y: &[f64]) -> f64 {

    let avg_x: f64 = avg(x);
    let avg_y: f64 = avg(y);
    let mut ret: f64 = 0.0;
    for (ix, iy) in zip(x, y) { ret += (ix - avg_x) * (iy - avg_y); }
    return ret / x.len() as f64;
}

/// Calculate weighted covariance from two columns of data of a column of weights
pub fn cov_weighted(x: &[f64], y: &[f64], weights: &[f64]) -> f64 {

    let avg_x: f64 = avg_weighted(x, weights);
    let avg_y: f64 = avg_weighted(y, weights);
    let mut cov: f64 = 0.0;
    for i in 0..weights.len() {
        cov += weights[i] * (x[i] - avg_x) * (y[i] - avg_y);
    }
    return cov / weights.iter().sum::<f64>();
}

/// Provides on-line statistics for N-dimensional samples
///
/// This struct accumulates observations without actually storing them, and on the fly provides
/// basic descriptive parameters of the accumulated sample.
///
/// See [`OnlineMultivariateWeighted`](OnlineMultivariateWeighted) for weighted statistics.
///
/// # Examples
/// The following example:
/// - creates a 2-dimensional normal distribution
/// - creates an `OnlineMultivariateStatistics` object
/// - withdraws 50000 observations from the distribution and accumulates them
/// - checks the convergence of expectations vector
///
/// Refer to [`MultiNormalDistribution`](bioshell_statistics::MultiNormalDistribution) struct documentation to see
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

/// Provides on-line statistics for N-dimensional weighted samples
///
/// This struct accumulates observations without actually storing them, and on the fly provides
/// basic descriptive parameters of the accumulated weighted sample.
///
/// See [`OnlineMultivariateStatistics`](OnlineMultivariateStatistics) for un-weighted statistics
///
/// # Examples
/// ```
/// use bioshell_statistics::{OnlineMultivariateWeighted, avg_weighted, var_weighted};
/// # use rand::rngs::SmallRng;
/// # use rand::{Rng, SeedableRng};
/// # use std::iter::zip;
/// // --- Here are 2D measurements and appropriate weights
/// let data_x = [3.0, 6.0, 5.0, 4.0, 7.0];
/// let data_y = [1.5, 3.0, 2.5, 2.0, 3.5];
/// let wghts = [0.2, 0.3, 0.4, 0.5, 0.6];
/// let mut stats = OnlineMultivariateWeighted::new(2);
/// // --- accumulate observations
/// for i in 0..data_x.len() {
///     stats.accumulate(&[data_x[i],data_y[i]], wghts[i]);
/// }
/// // --- check if statistics are correct
/// assert!((stats.avg()[0] - 5.3).abs() < 0.00001);
/// assert!((stats.avg()[1] - 2.65).abs() < 0.00001);
/// assert!((stats.var()[0] - 1.91).abs() < 0.00001);
///
/// // --- Generate thousand random 1D observations from [0,10) and respective weights from [0, 1.0)
/// let mut rng = SmallRng::seed_from_u64(0);
/// let sample: Vec<f64> = (0..100).map(|_| rng.gen_range(0.0..10.0)).collect();
/// let weights: Vec<f64> = (0..100).map(|_| rng.gen_range(0.0..1.0)).collect();
/// let mut stats = OnlineMultivariateWeighted::new(1);
/// for (s, w) in zip(&sample, &weights) { stats.accumulate(&[*s], *w);}
///
/// assert!((stats.avg()[0] - avg_weighted(&sample, &weights)).abs() < 0.001);
/// assert!((stats.var()[0] - var_weighted(&sample, &weights)).abs() < 0.001);
/// ```
#[derive(Clone)]
pub struct OnlineMultivariateWeighted {
    dim: usize,
    total: f64,
    m1: Vec<f64>,
    m2: Vec<f64>,
    min: Vec<f64>,
    max: Vec<f64>,
    cov: Vec<Vec<f64>>,
}

impl OnlineMultivariateWeighted {

    /// Create a new object to gather statistics on weighted N-dimensional samples
    ///
    /// # Arguments
    /// * `dim` - dimensionality of the data being observed
    ///
    pub fn new(dim: usize) -> OnlineMultivariateWeighted {
        OnlineMultivariateWeighted{dim, total:0.0, m1: vec![0.0; dim], m2: vec![0.0; dim],
            min: vec![0.0; dim], max: vec![0.0; dim],
            cov: vec![vec![0.0; dim]; dim]}
    }

    /// Returns the number of dimensions of the observed vectors
    pub fn dim(&self) -> usize { self.dim }

    /// Accumulate a single N-dimensional point
    pub fn accumulate(&mut self, d:&[f64], weight: f64) {

        assert_eq!(d.len(), self.dim);                  // --- incoming vector must be of the same size at the statistics
        if weight < 1.0e-10 { return;}
        if self.total < 0.000000000001 {                // --- if this is the very first observation...
            for i in 0..self.dim {
                self.min[i] = d[i];                     // --- copy it as min and max
                self.max[i] = d[i];
            }
        }
        self.total += weight;
        for i in 0..self.dim {
            let delta_x: f64 = d[i] - self.m1[i];
            self.m1[i] += delta_x * weight / self.total;  // --- M1[i] is now the new average for i-th dimension
            self.m2[i] += weight * delta_x * (d[i] -self.m1[i]);

            if self.m1[i].is_nan() {
                println!("{} {}   {} {:?} {:?}", &d[i], weight, self.total, &self.m1, &self.m2);
            }

            for j in i+1..self.dim {             // --- here we update the covariance terms
                let delta_y: f64 = d[j] - self.m1[j];
                self.cov[i][j] += weight * delta_y * (d[i] - self.m1[i]);
                self.cov[j][i] = self.cov[i][j];
            }

            self.min[i] = self.min[i].min(d[i]);
            self.max[i] = self.max[i].max(d[i]);
        }
    }

    /// Accumulate a single 1-dimensional point
    pub fn accumulate_1d(&mut self, x:f64, w: f64) {
        let v = [x];
        self.accumulate(&v, w);
    }

    /// Returns the sum of weight for the observed samples
    pub fn total(&self) -> f64{ self.total }

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
            ret[i] /= self.total;
        }

        return ret;
    }

    /// Returns the covariance matrix
    pub fn cov(&self) -> Vec<Vec<f64>> {

        let mut ret = self.cov.clone();
        for i in 0..ret.len() {
            for j in 0..i {
                ret[i][j] /= self.total;
                ret[j][i] /= self.total;
            }
            ret[i][i] = self.m2[i] / (self.total);
        }

        return ret;
    }
}

/// Calculates a given quantile of an arbitrary large sample of data.
///
/// The QuantileP2 struct uses the [,,P-Square Algorithm for Dynamic Calculation of Quantiles and Histograms without Storing Observations''](https://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf)
/// statistical algorithm to compute a desired quantile of data without actually storing them.
///
/// ```rust
///     let n_samples = 100000;
///
///     let normal = Normal::new(2.0, 3.0).unwrap();         // --- mean 2.0, standard deviation 3.0
///     let mut q2 = QuantileP2::new(0.5);
///     let mut q3 = QuantileP2::new(0.75);
///
///     for _ in 0..n_samples {
///         let x = normal.sample(&mut rand::thread_rng());
///         q2.accumulate(x);
///         q3.accumulate(x);
///     }
///     assert!((q2.quantile() - 2.0).abs() < 0.1);
///     assert!((q3.quantile() - 4.02347).abs() < 0.1);
///```
pub struct QuantileP2 {
    dn: [f64;5],        // Precalculated   marker increments
    npos: [f64;5],      // desired marker positions
    pos: [f64;5],       // Marker positions
    q: Vec<f64>,        // Marker heights
}

impl QuantileP2 {
    /// Start calculations of a desired p-quantile.
    ///
    /// At every [`accumulate()`](accumulate()) call, the object will improve it's estimate for the
    /// p-quantile, i.e. such a value `$x_p$` that `$P_X((-\infty, x_p]) \geqslant p$` and `$P_X([x_p,\infty)) \leqslant 1-p$`
    /// where `$P_X$` is a given (empirical) probability distribution.
    /// For example, `$x_{0.5}$` is the mean of the distribution.
    pub fn new(p: f64) -> QuantileP2 {
        return QuantileP2 {
            dn: [0.0, p/2.0, p, (1.0 + p)/2.0, 1.0],
            npos: [1.0, 1.0 + 2.0*p, 1.0 + 4.0*p, 3.0 + 2.0*p, 5.0],
            pos: [1.0, 2.0, 3.0, 4.0, 5.0],
            q: vec![],
        };
    }

    /// Accumulates the next point from a distribution and updates the current estimate of the quantile
    pub fn accumulate(&mut self, x: f64) {
        if self.q.len() <= 5 {
            self.q.push(x);
            self.q.sort_by(|a, b| a.partial_cmp(b).unwrap());
            if self.q.len() < 5 { return; }
        }

        let mut k: usize;
        if x < self.q[0] { // Update minimum value
            self.q[0] = x;
            k = 0;
        } else if self.q[0] <= x && x < self.q[1] { k = 0; }
        else if self.q[1] <= x && x < self.q[2] { k = 1; }
        else if self.q[2] <= x && x < self.q[3] { k = 2; }
        else if self.q[3] <= x && x <= self.q[4] { k = 3; }
        else { // Update maximum value
            self.q[4] = x;
            k = 3;
        }

        // Increment all positions starting at marker k+1
        for i in k + 1..5 { self.pos[i] += 1.0; }

        // Update desired marker positions
        for i in 0..5 { self.npos[i] += self.dn[i]; }

        // Adjust marker heights 2-4 if necessary
        for i in 1..4 {
            let d: f64 = self.npos[i] - self.pos[i];

            if (d >= 1.0 && (self.pos[i + 1] - self.pos[i]) > 1.0) || (d <= -1.0 && (self.pos[i - 1] - self.pos[i]) < -1.0) {
                let ds = d.signum();

                // Try adjusting q using P-squared formula
                let tmp: f64 = self.parabolic(ds, i);
                if self.q[i - 1] < tmp && tmp < self.q[i + 1] {
                    self.q[i] = tmp;
                } else {
                    self.q[i] = self.linear(ds, i);
                }
                self.pos[i] += ds;
            }
        }
    }

    /// Returns the current estimate of the quantile
    pub fn quantile(&self) -> f64 { self.q[2] }

    fn linear(&self, d: f64, i: usize) -> f64 {
        return self.q[i] + d * (self.q[i + d as usize] - self.q[i]) / (self.pos[i + d as usize] - self.pos[i]);
    }

    fn parabolic(&self, d: f64, i: usize) -> f64 {
        let a: f64 = d / (self.pos[i + 1] - self.pos[i - 1]);

        let b: f64 = (self.pos[i] - self.pos[i - 1] + d) * (self.q[i + 1] - self.q[i]) / (self.pos[i + 1] - self.pos[i])
            + (self.pos[i + 1] - self.pos[i] - d) * (self.q[i] - self.q[i - 1]) / (self.pos[i] - self.pos[i - 1]);

        return self.q[i] + a * b;
    }
}

