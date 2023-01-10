mod histograms;
mod distributions;
mod descriptive;

pub use descriptive::{OnlineMultivariateStatistics};

pub use histograms::{Histogram};
pub use distributions::{NormalDistribution, MultiNormalDistribution, Estimable, Distribution};