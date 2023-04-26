mod histograms;
mod distributions;
mod descriptive;

pub use descriptive::{OnlineMultivariateStatistics, OnlineMultivariateWeighted,
                      avg, avg_weighted, var, var_weighted, cov, cov_weighted, avg_var_weighted};

pub use histograms::{Histogram};
pub use distributions::{NormalDistribution, MultiNormalDistribution, Estimable, Distribution};