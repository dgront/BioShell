use std::fmt;
use nalgebra::{DMatrix, DVector};
use rand::Rng;
use rand_distr::Distribution as OtherDistribution;

use crate::{OnlineMultivariateStatistics};

// ========== Distribution trait ==========

/// Defines what any probability distribution must offer.
///
/// A probability distribution function (pdf) trait primarily requires that a derived type will be able to:
///   - evaluate the pdf(x) at any given ``x``
///   - draw a random sample from the distribution
/// Since this trait covers also multivariate distributions, the ``x`` is a vector, i.e. ``Vec<f64>``
pub trait Distribution {
    /// Evaluates the probability distribution function at a given point
    fn pdf(&self, x: &Vec<f64>) -> f64;

    /// Withdraws a random observation from this probability distribution
    fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R, out: &mut [f64]);

    /// Says how many dimensions this distribution has
    fn dim(&self) -> usize;
}


// ========== Normal probability distribution structure and implementation ==========
/// Normal probability distribution
///
/// ```
/// use bioshell_statistics::{Distribution, NormalDistribution};
/// // --- create a N(2.0, 1.5) distribution
/// let mut n: NormalDistribution = NormalDistribution::new(2.0, 1.5);
///
/// // --- evaluate pdf at x = 1.0
/// let prob = n.pdf(&vec![1.0]);
/// assert!((prob - 0.21297).abs() < 0.0001);   // 0.21297 is the value from statistical tables
/// ```
#[derive(Clone)]
pub struct NormalDistribution {
    mean: f64,
    sigma: f64,
    const_1: f64,
    normal_generator: rand_distr::Normal<f64>            // for rnd sampling
}

// #[pymethods]
impl NormalDistribution {
    /// Creates a new normal probability distribution N(mu, sigma)
    // #[new]
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

    fn sample<R: Rng + ?Sized>(&mut self, rand: &mut R, out: &mut [f64]) {
        out[0] = self.normal_generator.sample(rand);
    }

    fn dim(&self) -> usize { return 1;  }
}

// ========== Multivariate Normal probability distribution structure and implementation ==========

/// N-dimentional normal probability distribution
///
/// # Examples
/// To properly initialize an N-dimensional Gaussian distribution, one has to prepare a vector
/// of expected values and a covariance matrix; both are objects of types defined in the
/// `nalgebra` crate: `DMatrix` and `DVector`
/// ```
/// use bioshell_statistics::{Distribution, MultiNormalDistribution};
/// use nalgebra::{DMatrix, DVector};
///
/// // --- create a 2-D normal distribution
/// let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
/// // --- set expected values and a covariance matrix
/// n.set_parameters(&DVector::from_vec(vec![0.1, 0.1]),
///         &DMatrix::from_vec(2, 2, vec![0.2, 0.1, 0.1, 2.0]));
///
/// // --- evaluate pdf at x = [1.0, 0.0]
/// let logprob = n.logpdf(&vec![1.0, 0.0]);
/// assert!((logprob + 3.469636899044226).abs() < 0.0001);
/// ```
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

    /// Returns log-probability for a given pdf vector
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
            self.mean[i] = stats.avg()[i];
            let covar = stats.cov();
            for j in 0..self.dim {
                self.sigma[(i, j)] = covar[i][j];
            }
            self.sigma[(i, i)] = stats.var()[i];
        }
        self.setup();
    }
}

impl Distribution for MultiNormalDistribution {

    fn pdf(&self, x: &Vec<f64>) -> f64 {
        self.logpdf(x).exp()
    }

    fn sample<R: Rng + ?Sized>(&mut self, rng: &mut R, out: &mut [f64]) {
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

/// Provides an ability to estimate parameters of a probability [`Distribution`](Distribution)
///
/// Any [`Distribution`](Distribution) that is [`Estimable`](Estimable), must implement
/// [`estimate()`](estimate) method, which estimates parameters of this distribution from
/// a provided data sample. Resulting parameters are stored within `self` i.e. the
///  `Distribution` object estimates itself from the sample
pub trait Estimable {
    /// Estimate parameters of a probability distribution from a given sample
    ///
    /// # Arguments
    /// * `sample` - observations used to infer the distribution parameters;
    /// dimensionality of the data must match this `Distribution` object
    fn estimate(&mut self, sample: &Vec<Vec<f64>>);

    /// Estimate parameters of a probability distribution from a given sample fraction
    ///
    /// # Arguments
    /// * `sample` - observations used to infer the distribution parameters
    /// * `selection` - boolean flags pointing to which data rows of a given sample should actually be used
    ///     for estimation
    fn estimate_from_selected(&mut self, sample: &Vec<Vec<f64>>, selection: &Vec<bool>);
}

impl Estimable for NormalDistribution {

    /// Estimate parameters of this [`NormalDistribution`](NormalDistribution) from a given sample
    ///
    /// # Example
    /// ```
    /// use bioshell_statistics::{Estimable, NormalDistribution};
    ///
    /// let mut n: NormalDistribution = NormalDistribution::new(0.0, 1.0);
    ///
    /// // --- evaluate pdf at x = 1.0
    /// let sample = vec![vec![1.1], vec![1.2], vec![1.4], vec![1.3]];
    /// n.estimate(&sample);
    /// assert!((n.mean() - 1.25).abs() < 0.0001);
    /// ```
    fn estimate(&mut self, sample: &Vec<Vec<f64>>) {
        assert_eq!(sample[0].len(), 1);

        let mut stats = OnlineMultivariateStatistics::new(1);
        sample.iter().for_each(|x| stats.accumulate(x));
        self.set_parameters(stats.avg()[0], stats.var()[0].sqrt());
    }

    fn estimate_from_selected(&mut self, sample: &Vec<Vec<f64>>, selection: &Vec<bool>) {
        assert_eq!(sample[0].len(), 1);

        let mut stats = OnlineMultivariateStatistics::new(1);
        for i in 0..sample.len() {
            if selection[i] { stats.accumulate(&sample[i]);  }
        }
        self.set_parameters(stats.avg()[0], stats.var()[0].sqrt());
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
    /// use std::fmt::Write;
    /// use bioshell_statistics::NormalDistribution;
    ///
    /// let nd = NormalDistribution::new(2.0, 0.5);
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of the distribution
    /// write!(actual, "{}", nd).unwrap();
    ///
    /// let expected = "N(mu= 2.0000, sdev= 0.5000)";
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
    ///
    /// # Examples
    /// Create a `MultiNormalDistribution` and turn it into a string
    /// ```rust
    /// use std::fmt::Write;
    /// use nalgebra::{DMatrix, DVector};
    /// use bioshell_statistics::MultiNormalDistribution;
    ///
    /// let mut n: MultiNormalDistribution = MultiNormalDistribution::new(2);
    /// n.set_parameters(&DVector::from_vec(vec![1.0, 2.0]),
    /// &DMatrix::from_vec(2,2, vec![1.0, 0.5, 0.5, 1.0]));
    ///
    /// let expected = "'mu': [ 1.0000,  2.0000], 'sigma': [ [ 1.0000,  0.5000], [ 0.5000,  1.0000]]";
    /// let mut actual = String::new();
    /// write!(actual, "{}", n).unwrap();
    /// assert_eq!(actual, expected);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::from("'mu':");
        print_vector!(&self.mean(), out);
        // print_vector!(&self.sigma()[i], out);
        out.push_str(", 'sigma': [");
        print_vector!(&self.sigma().row(0), out);
        for i in 1..self.dim() {
            out.push_str(",");
            print_vector!(&self.sigma().row(i), out);
        }
        out.push_str("]");

        write!(f, "{}", out)
    }
}


