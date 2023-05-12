//! Staple statistical utilities necessary in daily life.
//!
//! This crate provides basic functions to compute descriptive statistics, such as [`avg()`](avg())
//! or [`var_weighted()`](var_weighted()), which accept a vector of observations as an argument.
//!
//!
//! Such statistics can also be evaluated using an on-line calculator:
//! [`OnlineMultivariateStatistics`](OnlineMultivariateStatistics) in the regular case
//! and [`OnlineMultivariateWeighted`](OnlineMultivariateWeighted) for weighted sample.
//! The structs solve the problem of generating accurate statistics when a data set is too large to fit in memory.
//! One can also create a [`Histogram`](Histogram) of arbitrary data sample.
//!
//! Finally, the crate provides [`Distribution`](Distribution) trait; a few basic distributions have been already implemented.
//! In particular,  [`Estimable`](Estimable) distributions are used by expectation-maximization
//! algorithm implemented by ``bioshell-clustering`` crate
//!
mod histograms;
mod distributions;
mod descriptive;

pub use descriptive::{OnlineMultivariateStatistics, OnlineMultivariateWeighted, QuantileP2,
                      avg, avg_weighted, var, var_weighted, cov, cov_weighted, avg_var_weighted};

pub use histograms::{Histogram};
pub use distributions::{NormalDistribution, MultiNormalDistribution, Estimable, Distribution};