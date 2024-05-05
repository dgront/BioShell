//! Provides distance functions
//!
//! A distance function computes a distance between two points of generic type `T`,
//! which must be indexable, i.e. must be bounded by `Index<usize, Output = f64>` trait.
//! The signature of each distance function is `Fn(&T, &T, usize) -> f64` i.e. it accepts
//! two arguments of type `&T` and `usize` the number of dimensions. The latter tells
//! what is the maximum value of an index used to reach coordinates of `T` points

use std::ops::Index;

/// Calculate the squared euclidean distance between two points.
///
/// The arguments of a generic type `T` must provide the indexing operator with returns `f64`
/// values of a coordinate of that point.
///
/// # Arguments
/// * `a` - the first k-dimensional point
/// * `b` - the second k-dimensional point
/// * `dimensionality` - the number of dimensions for each of the two points
///
/// ```rust
/// use bioshell_datastructures::euclidean_distance_squared;
/// let d = euclidean_distance_squared(&[0.1, 0.1], &[0.2, 0.2], 2);
/// assert!((d-0.02).abs() < 0.000001);
/// let d = euclidean_distance_squared(&vec![0.1, 0.1], &vec![0.2, 0.2], 2);
/// assert!((d-0.02).abs() < 0.000001);
/// ```
pub fn euclidean_distance_squared<T>(a: &T, b: &T, dimensionality: usize) -> f64 where T: Index<usize, Output = f64> {

    let mut ret = 0.0;
    for i in 0..dimensionality { ret += (a[i] - b[i]) * (a[i] - b[i]); }
    return ret;
}

