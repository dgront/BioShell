use std::fmt::Display;

// use pyo3::prelude::*;

use bioshell_statistics::{Distribution, Estimable};

// ========== Expectation - Maximization ==========
pub fn expectation_maximization<D: Estimable + Distribution + Display>(
    distributions: &mut Vec<D>,
    data: &Vec<Vec<f64>>,
    assignment: &mut Vec<Vec<bool>>,
    epsilon: f64,
) -> f64 {
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
        progress = ((log_liklhd - last_log_liklhd) / log_liklhd).abs();
        last_log_liklhd = log_liklhd;
        // println!("{} {}",last_log_liklhd, progress);
    }

    return last_log_liklhd;
}

fn which_distribution<D: Estimable + Distribution>(
    distributions: &Vec<D>,
    p: &Vec<f64>,
) -> (usize, f64) {
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

// #[pymodule]
// fn bioshell_numerical(_: Python, m: &PyModule) -> PyResult<()> {
//     m.add_class::<NormalDistribution>()?;
//     Ok(())
// }
