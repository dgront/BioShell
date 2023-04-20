//!Provides distance functions

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
/// use bioshell_numerical::distance::euclidean_distance_squared;
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