/// Provides on-line statistics for  N-dimensional samples
pub struct OnlineMultivariateStatistics {
    dim: usize,
    count: usize,
    m1: Vec<f64>,
    m2: Vec<f64>,
    cov: Vec<Vec<f64>>,
}

impl OnlineMultivariateStatistics {

    /// Create a new object to gather statistics on N-dimensional samples
    pub fn new(dim: usize) -> OnlineMultivariateStatistics {
        OnlineMultivariateStatistics{dim, count:0, m1: vec![0.0; dim], m2: vec![0.0; dim],
            cov: vec![vec![0.0; dim]; dim]}
    }

    /// Returns the dimension od the observed vectors
    pub fn dim(&self) -> usize { self.dim }

    /// Accumulate a single N-dimensional point
    pub fn accumulate(&mut self, d:&Vec<f64>) {

        assert_eq!(d.len(), self.dim);                  // --- incoming vector must be of the same size at the statistics

        self.count += 1;
        for i in 0..self.dim {
            let delta_x: f64 = d[i] - self.m1[i];
            self.m1[i] += delta_x / self.count as f64;  // --- M1[i] is now the new average for i-th dimension
            self.m2[i] += delta_x * (d[i] -self.m1[i]);

            for j in i+1..self.dim {
                let delta_y: f64 = d[j] - self.m1[j];
                self.cov[i][j] += delta_y * (d[i] - self.m1[i]);
                self.cov[j][i] = self.cov[i][j];
            }
        }
    }

    /// Returns the number of observed samples
    pub fn count(&self) ->usize{ self.count }

    /// Returns the average of the values observed so far.
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate
    pub fn avg(&self, id:usize) ->f64 { self.m1[id] }

    /// Returns the variance of the values observed so far.
    ///
    ///  # Arguments
    /// * `id` - index of the coordinate
    pub fn var(&self, id:usize) ->f64 { self.m2[id] / (self.count as f64 - 1.0) }

    /// Returns the covariance between i-th and j-th columns of the data observed so far
    ///
    ///  # Arguments
    /// * `i` - index of the coordinate
    /// * `j` - index of the coordinate
    pub fn covar(&self, i:usize, j:usize) ->f64 { self.cov[i][j] / (self.count as f64 - 1.0) }
}