use std::fmt;
use std::fmt::Display;
use nalgebra::{DMatrix, DVector};
use rand::Rng;
use rand_distr::Distribution as OtherDistribution;

use crate::statistics::OnlineMultivariateStatistics;

// ========== Distribution trait ==========

/// Defines what any probability distribution must offer
pub trait Distribution {
    /// Evaluates the probability distribution function at a given point
    fn pdf(&self, x: &Vec<f64>) -> f64;

    /// Withdraws a random observation from this probability distribution
    fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R, out: &mut Vec<f64>);

    /// Says how many dimensions this distribution has
    fn dim(&self) -> usize;
}


// ========== Normal probability distribution structure and implementation ==========
/// Normal probability distribution
#[derive(Clone)]
pub struct NormalDistribution {
    mean: f64,
    sigma: f64,
    const_1: f64,
    normal_generator: rand_distr::Normal<f64>            // for rnd sampling
}

impl NormalDistribution {
    /// Creates a new normal probability distribution N(mu, sigma)
    pub fn new(mu: f64, sigma: f64) -> NormalDistribution {
        NormalDistribution { mean: mu, sigma,
            const_1: (1.0 / (sigma * (std::f64::consts::PI * 2.0).sqrt())),
            normal_generator: rand_distr::Normal::new(mu, sigma).unwrap()
        }
    }

    /// Returns the mean (expected) value of this distribution
    pub fn mean(&self) -> f64 { self.mean }

    /// Returns the standard deviation of this distribution
    pub fn sdev(&self) -> f64 { self.sigma }

    /// Sets parameters of this normal probability distribution to mu and sigma
    pub fn set_parameters(&mut self, mu: f64, sigma: f64) {
        self.mean = mu;
        self.sigma = sigma;
        self.const_1 = 1.0 / (sigma * (std::f64::consts::PI * 2.0).sqrt());
        self.normal_generator = rand_distr::Normal::new(mu, sigma).unwrap();
    }
}

impl Distribution for NormalDistribution {
    fn pdf(&self, x: &Vec<f64>) -> f64 {
        let ix = x[0] - self.mean;
        return self.const_1 * (-(ix * ix) / (2.0 * self.sigma * self.sigma)).exp();
    }

    fn sample<R: Rng + ?Sized>(&mut self, rand: &mut R, out: &mut Vec<f64>) {
        out[0] = self.normal_generator.sample(rand);
    }

    fn dim(&self) -> usize { return 1;  }
}

// ========== Multivariate Normal probability distribution structure and implementation ==========

#[derive(Clone)]
pub struct MultiNormalDistribution {
    dim: usize,
    logdet: f64,
    mean: DVector<f64>,
    sigma: DMatrix<f64>,
    u: DMatrix<f64>,
    cku: DMatrix<f64>,                                  // Cholesky decomposition to random sampling from the distribution
    tmp: DVector<f64>,                                  // for rnd sampling
    normal_generator: rand_distr::Normal<f64>           // for rnd sampling
}

impl MultiNormalDistribution {

    /// Creates a new multivariate normal distribution
    /// Parameters are set to a vector of zeros and a unity matrix (expected and covariance, respectively)
    pub fn new(dim: usize) -> MultiNormalDistribution {

        let mu = DVector::<f64>::zeros(dim);
        let sig = DMatrix::<f64>::identity(dim, dim);
        let u = DMatrix::<f64>::identity(dim, dim);
        let cku = DMatrix::<f64>::identity(dim, dim);
        let tmp = DVector::<f64>::zeros(dim);
        let mut out = MultiNormalDistribution { dim, logdet:0.0 , mean: mu, sigma: sig, u, cku, tmp,
            normal_generator: rand_distr::Normal::new(0.0, 1.0).unwrap()};
        out.setup();

        return out;
    }

    /// Returns the immutable reference to the mean (expected) vector of this distribution
    pub fn mean(&self) -> &DVector<f64> { &self.mean }

    /// Returns the immutable reference to the standard deviation matrix of this distribution
    pub fn sigma(&self) -> &DMatrix<f64> { &self.sigma }

    /// Sets parameters of this multinomial normal distribution
    pub fn set_parameters(&mut self, mu:&DVector::<f64>, sigma: &DMatrix::<f64>) {
        self.mean = mu.clone();
        self.sigma = sigma.clone();
        self.setup();
    }

    pub fn logpdf(&self, x: &Vec<f64>) -> f64 {
        let n: usize = x.len();
        let mut xm = DVector::from_vec(x.clone());
        xm -= &self.mean;
        let mut dp = DVector::<f64>::zeros(n);
        for i in 0..n {
            dp[i] = xm.dot(&self.u.column(i));
        }
        let log2pi     = (2.0*std::f64::consts::PI).ln();
        let maha_dist: f64 = dp.iter().map(|x| x*x).sum();
        let val = -0.5 * (self.dim as f64 * log2pi + maha_dist + self.logdet);

        return val;
    }


    /// Sets up internal data ready for PDF evaluation
    /// see https://gregorygundersen.com/blog/2019/10/30/scipy-multivariate/
    fn setup(&mut self) {
        let sig_copy = self.sigma.clone();
        let eigen = sig_copy.symmetric_eigen();

        self.logdet = eigen.eigenvalues.map(|e| e.ln()).sum();
        let mut valsinv: DVector<f64> = DVector::<f64>::zeros(self.dim);
        let mut valsinv_sqrt: DVector<f64> = DVector::<f64>::zeros(self.dim);
        for i in 0..self.dim {
            valsinv[i] = 1.0 / eigen.eigenvalues[i];
            valsinv_sqrt[i] = valsinv[i].sqrt();
        }
        self.u = eigen.eigenvectors.clone();
        for i in 0..self.dim {
            for j in 0..self.dim {
                self.u[(i,j)] = self.u[(i,j)] * valsinv_sqrt[j];
            }
        }

        self.cku = self.sigma.clone().cholesky().unwrap().l();
    }

    /// Sets mean and sigma parameters from observations gathered
    ///in a given  OnlineMultivariateStatistics object
    fn set_from_statistics(&mut self, stats: &OnlineMultivariateStatistics) {

        assert_eq!(stats.dim(), self.dim);

        for i in 0..self.dim {
            self.mean[i] = stats.avg(i);
            for j in 0..self.dim {
                self.sigma[(i, j)] = stats.covar(i, j);
            }
            self.sigma[(i, i)] = stats.var(i);
        }
        self.setup();
    }
}

impl Distribution for MultiNormalDistribution {

    fn pdf(&self, x: &Vec<f64>) -> f64 {
        self.logpdf(x).exp()
    }

    fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R, out: &mut Vec<f64>) {
        for i in 0..self.dim {
            self.tmp[i] = self.normal_generator.sample(rng);
        }
        let o = &self.cku * &self.tmp + &self.mean;
        for i in 0..self.dim {
            out[i] = o[i];
        }
    }

    fn dim(&self) -> usize { return self.dim;  }
}

// ========== Estimable trait and its implementation for distributions ==========
pub trait Estimable {
    fn estimate(&mut self, sample: &Vec<Vec<f64>>);
    fn estimate_from_selected(&mut self, sample: &Vec<Vec<f64>>, selection: &Vec<bool>);
}

impl Estimable for NormalDistribution {

    fn estimate(&mut self, sample: &Vec<Vec<f64>>) {
        assert_eq!(sample[0].len(), 1);

        let mut stats = OnlineMultivariateStatistics::new(1);
        sample.iter().for_each(|x| stats.accumulate(x));
        self.set_parameters(stats.avg(0), stats.var(0).sqrt());
    }

    fn estimate_from_selected(&mut self, sample: &Vec<Vec<f64>>, selection: &Vec<bool>) {
        assert_eq!(sample[0].len(), 1);

        let mut stats = OnlineMultivariateStatistics::new(1);
        for i in 0..sample.len() {
            if selection[i] { stats.accumulate(&sample[i]);  }
        }
        self.set_parameters(stats.avg(0), stats.var(0).sqrt());
    }
}

impl Estimable for MultiNormalDistribution {

    fn estimate(&mut self, sample: &Vec<Vec<f64>>) {
        assert_eq!(sample[0].len(), self.dim);

        let mut stats = OnlineMultivariateStatistics::new(self.dim);
        sample.iter().for_each(|x| stats.accumulate(x));

        self.set_from_statistics(&stats);
    }

    fn estimate_from_selected(&mut self, sample: &Vec<Vec<f64>>, selection: &Vec<bool>) {
        assert_eq!(sample[0].len(), self.dim);

        let mut stats = OnlineMultivariateStatistics::new(self.dim);
        for i in 0..sample.len() {
            if selection[i] { stats.accumulate(&sample[i]);  }
        }

        self.set_from_statistics(&stats);
    }
}

// ========== Printing distributions nicely ==========

impl fmt::Display for NormalDistribution {
    /// Nicely prints a NormalDistribution object
    /// # Examples
    ///
    /// Create a `NormalDistribution` and turn it into a string
    ///
    /// ```rust
    /// use bioshell_core::Sequence;
    /// use std::fmt::Write;
    ///
    /// let nd = NormalDistribution::new(2.0, 0.5);
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of the distribution
    /// write!(actual, "{}", nd).unwrap();
    ///
    /// let expected = "N(mu=    5, sdev=  0.5)";
    /// assert_eq!(actual, expected)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "N(mu={:7.4}, sdev={:7.4})", self.mean(), self.sdev())
    }
}

macro_rules!  print_vector {
    ($vec:expr, $out_str:expr) => {
        $out_str = format!("{} [{:7.4}", $out_str, $vec[0]);
        for i in 1..$vec.len() { $out_str = format!("{}, {:7.4}", $out_str, $vec[i]);}
        $out_str.push_str("]");
    }
}

impl fmt::Display for MultiNormalDistribution {
    /// Nicely prints a MultiNormalDistribution object
    /// # Examples
    ///
    /// Create a `MultiNormalDistribution` and turn it into a string
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::from("mu = ");
        print_vector!(&self.mean(), out);
        // print_vector!(&self.sigma()[i], out);
        out.push_str(", sigma = [");
        print_vector!(&self.sigma().row(0), out);
        for i in 1..self.dim() {
            out.push_str(",");
            print_vector!(&self.sigma().row(i), out);
        }
        out.push_str("]");

        write!(f, "{}", out)
    }
}

// ========== Expectation - Maximization ==========
pub fn expectation_maximization<D: Estimable+Distribution+Display>(distributions:&mut Vec<D>,
        data: &Vec<Vec<f64>>, assignment: &mut Vec<Vec<bool>>, epsilon: f64) -> f64 {

    // --- check whether the input vectors have the right size
    assert_eq!(data.len(), assignment[0].len());
    assert_eq!(distributions.len(), assignment.len());

    let n_data: usize = data.len();
    let n_dist: usize = distributions.len();
    let mut progress: f64 = 1.0;
    let mut last_log_liklhd = 0.0;

    while progress > epsilon {
        // ---  The "expectation" step: assignment of points to the best distributions
        for i in 0..n_dist {
            distributions[i].estimate_from_selected(&data, &assignment[i]);
            assignment[i].fill(false);
        }
        // ---  The "maximization" step: assignment of points to the best distributions
        let mut log_liklhd: f64 = 0.0;
        for i in 0..n_data {
            let (best_idx, best_val) = which_distribution(&distributions, &data[i]);
            assignment[best_idx][i] = true;
            log_liklhd += best_val.ln();
        }
        progress = ((log_liklhd - last_log_liklhd)/log_liklhd).abs();
        last_log_liklhd = log_liklhd;
        // println!("{} {}",last_log_liklhd, progress);
    }

    return last_log_liklhd;
}

fn which_distribution<D: Estimable+Distribution>(distributions: &Vec<D>, p: &Vec<f64>)  -> (usize, f64) {

    let mut best_idx: usize = 0;
    let mut best_val: f64 = distributions[0].pdf(p);
    for i in 1..distributions.len() {
        let v = distributions[i].pdf(p);
        if v > best_val {
            best_idx = i;
            best_val = v
        }
    }
    return (best_idx, best_val);
}