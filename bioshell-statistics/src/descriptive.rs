/// Provides on-line statistics for N-dimensional samples
///
/// This struct accumulates observations without actually storing them, and on the fly provides
/// basic descriptive parameters of the accumulated sample
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
    pub fn accumulate(&mut self, d:&Vec<f64>) {

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
    pub fn accumulate_1D(&mut self, x:f64) {
        let v = vec![x];
        self.accumulate(&v);
    }

    /// Returns the number of observed samples
    pub fn count(&self) -> usize{ self.count }

    /// Returns the minimum value observed in a given dimension
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate i.e. dimension (starts form 0)
    pub fn min(&self, id:usize) -> f64 { self.min[id] }

    /// Returns the maximum value observed in a given dimension
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate i.e. dimension (starts form 0)
    pub fn max(&self, id:usize) -> f64 { self.max[id] }

    /// Returns the average of the values observed so far.
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate i.e. dimension (starts form 0)
    pub fn avg(&self, id:usize) ->f64 { self.m1[id] }

    /// Returns the variance of the values observed so far.
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate i.e. dimension (starts form 0)
    pub fn var(&self, id:usize) ->f64 { self.m2[id] / (self.count as f64 - 1.0) }

    /// Returns the covariance between i-th and j-th columns of the data observed so far
    ///
    ///  # Arguments
    /// * `i` - index of the coordinate i.e. dimension (starts form 0)
    /// * `j` - index of the coordinate i.e. dimension (starts form 0)
    pub fn covar(&self, i:usize, j:usize) ->f64 { self.cov[i][j] / (self.count as f64 - 1.0) }
}