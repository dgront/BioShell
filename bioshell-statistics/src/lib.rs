mod descriptive;
mod distributions;
mod histograms;

pub use descriptive::OnlineMultivariateStatistics;

pub use distributions::{Distribution, Estimable, MultiNormalDistribution, NormalDistribution};
pub use histograms::Histogram;
