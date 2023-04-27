//! Expectation-Maximization algorithm fits a mixture of  generic distributions
//!
//!
use std::fmt::Display;
use std::iter::zip;
use log::{debug};

use bioshell_statistics::{Estimable, Distribution};


pub fn expectation_maximization<D: Estimable+Distribution+Display>(distributions:&mut Vec<D>,
        data: &Vec<Vec<f64>>, weights: &mut Vec<f64>, epsilon: f64) -> f64 {

    // --- check whether the input vectors have the right size
    assert_eq!(distributions.len(), weights.len());

    let n_data: usize = data.len();
    let n_dist: usize = distributions.len();
    let mut progress: f64 = 1.0;
    let mut last_log_liklhd = 0.0;

    let mut istep = 0;
    let mut reponsibity = vec![vec![0.0; n_data]; n_dist];
    while progress > epsilon {
        istep += 1;
        // ---  The "expectation" step: calculate responsibilities i.e. probability that i-th point belong to k-th distribution
        for k in 0..n_dist { weights[k] = 0.0}
        for i in 0..n_data {
            let p = &data[i];
            let mut total_r = 0.0;
            for k in 0..n_dist {
                let r= distributions[k].pdf(p);
                total_r += r;
                reponsibity[k][i] = r;
                weights[k] += r;
            }
            // --- normalize the responsibilities for each point (so the sum over distributions = 1.0)
            for k in 0..n_dist { reponsibity[k][i] /= total_r;}
        }
        let total_w: f64 = weights.iter().sum();
        // --- normalize the mixing weights
        for k in 0..n_dist { weights[k] /= total_w;}

        // ---  The "maximization" step: optimize parameters of each distribution
        for k in 0..n_dist {
            distributions[k].estimate_weighted(data, &reponsibity[k]);
        }
        // ---  Calculate log-likelihood to check the convergence
        let log_liklhd: f64 = log_likelihood(&distributions, &weights, data);
        progress = ((log_liklhd - last_log_liklhd) / log_liklhd).abs();
        last_log_liklhd = log_liklhd;
        debug!("log-likelihood after {} iteration: {}; change: {:.4}", istep, last_log_liklhd, progress);
    }

    return last_log_liklhd;
}

/// Calculates log-likelihood of fitness the given data to a mixture of probability distributions.
///
/// This function returns the log-likelihood value defined as:
/// ```math
/// \mathcal{L} = \sum_{i=1}^{N} \log \bigl( \sum_{c=1}^{K} w_c pdf_c(x_i) \bigr)
/// ```
/// which describes how well a set of `N` observed data points fits to a mixture of `K` distributions.
/// `$w_c$` denotes a weight of `c` component (distributions), `$1 \le c \le K$`.
pub fn log_likelihood<D: Estimable+Distribution+Display>(distributions:&Vec<D>,
         weights: &Vec<f64>, data: &Vec<Vec<f64>>) -> f64 {

    let mut log_liklhd: f64 = 0.0;
    for i in 0..data.len() {
        let mut point_liklhd = 0.0;
        // zip(weights, distributions).map(|w, d| point_liklhd += w* d.pdf(&data[i]));
        for (w, d) in zip(weights, distributions) {
            point_liklhd += w * d.pdf(&data[i]);
        }
        log_liklhd += point_liklhd.ln();
    }
    return log_liklhd;
}