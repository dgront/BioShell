use nalgebra::{DMatrix, DVector};
use rand_distr::{Normal};
use rand_distr::Distribution as OtherDistribution;

use crate::statistics::OnlineMultivariateStatistics;

// ========== Distribution trait ==========

pub trait Distribution {
    fn pdf(&self, x: &Vec<f64>) -> f64;
}


// ========== Normal probability distribution structure and implementation ==========
pub struct NormalDistribution {
    mean: f64,
    sigma: f64,
    const_1: f64
}

impl NormalDistribution {
    pub fn new(mu: f64, sigma: f64) -> NormalDistribution {
        NormalDistribution { mean: mu, sigma, const_1: (1.0 / (sigma * (std::f64::consts::PI * 2.0).sqrt())) }
    }
}

impl Distribution for NormalDistribution {
    fn pdf(&self, x: &Vec<f64>) -> f64 {
        let ix = x[0] - self.mean;
        return self.const_1 * (-(ix * ix) / (2.0 * self.sigma * self.sigma)).exp();
    }
}

// ========== Multivariate Normal probability distribution structure and implementation ==========

pub struct MultiNormalDistribution {
    dim: usize,
    logdet: f64,
    mean: DVector<f64>,
    sigma: DMatrix<f64>,
    u: DMatrix<f64>,
    cku: DMatrix<f64>,                                  // Cholesky decomposition to random sampling from the distribution
    tmp: DVector<f64>,                                  // for rnd sampling
    normal_generator:rand_distr::Normal<f64>            // for rnd sampling
}

impl MultiNormalDistribution {
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

    pub fn rand(&mut self, out: &mut Vec<f64>) {
        for i in 0..self.dim {
            self.tmp[i] = self.normal_generator.sample(&mut rand::thread_rng());
        }
        let o = &self.cku * &self.tmp + &self.mean;
        for i in 0..self.dim {
            out[i] = o[i];
        }
    }

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
}

impl Distribution for MultiNormalDistribution {

    fn pdf(&self, x: &Vec<f64>) -> f64 {
        self.logpdf(x).exp()
    }
}

// ========== Estimable trait and its implementation for distributions ==========
pub trait Estimable {
    fn estimate(&mut self, sample: &Vec<Vec<f64>>);
}


impl Estimable for MultiNormalDistribution {

    fn estimate(&mut self, sample: &Vec<Vec<f64>>) {
        assert_eq!(sample[0].len(), self.dim);

        let mut stats = OnlineMultivariateStatistics::new(self.dim);
        sample.iter().for_each(|x| stats.accumulate(x));

        for i in 0..self.dim{
            self.mean[i] = stats.avg(i);
            for j in 0..self.dim {
                self.sigma[(i,j)] = stats.covar(i,j);
            }
        }
        self.setup();
    }
}

// ========== Expectation - Maximization ==========
