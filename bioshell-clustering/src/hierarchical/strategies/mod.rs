//! Provides merging rules for hierarchical clustering.
//!
//! The merging rules are used to determine the distance between two clusters, based on the
//! distances between the clusters elements themselves. There are four merging rules implemented:
//!
//!  - `Single link`: the distance between two clusters is the **minimum** distance between any two elements of the two clusters.
//!  - `Complete link`: the distance between two clusters is the **maximum** distance between any two elements of the two clusters.
//!  - `Average link`: the distance between two clusters is the **average** distance between any two elements of the two clusters.
//!  - `Median link`: the distance between two clusters is based on the **median** link criterion
//!  - `Ward link`: the distance between two clusters is based on the **minimum variance** criterion (WPGMA)
//!  - `Centroid link`: the distance between two clusters is based on the centroid criterion
//!
//! These merging rules can be evaluated in constant time, which requires some additional information.
//! Specifically, the distance between any cluster `L` and a cluster `K`  just created by merging
//! two clusters `I` and `J`,  requires the sizes of the clusters `I`, `J` and `K`,
//! and the distances between the clusters `I` and `J`, `I` and `K`, and `J` and `K`.
//!

/// Calculates the distance between two clusters according to the `Single link` merging rule.
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = min(d(C_i,C_k),d(C_j,C_k))
/// ```
pub fn single_link(_size_i: usize, _size_j: usize, _size_k: usize,
                   _distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {
    distance_ik.min(distance_jk)
}

/// Calculates the distance between two clusters according to the `Complete link` merging rule.
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = max(d(C_i,C_k),d(C_j,C_k))
/// ```
pub fn complete_link(_size_i: usize, _size_j: usize, _size_k: usize,
                     _distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {
    distance_ik.max(distance_jk)
}

/// Calculates the distance between two clusters according to the `Complete link` merging rule.
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = \frac{|C_i|}{|C_i|+|C_j|}d(C_i,C_k) +  \frac{|C_j|}{|C_i|+|C_j|}d(C_j,C_k)
/// ```
pub fn average_link(size_i: usize, size_j: usize, _size_k: usize,
                    _distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {

    let d = 1.0 / ((size_i + size_j) as f32);

    return d * size_i as f32 * distance_ik + d * size_j as f32  * distance_jk;
}

/// Calculates the distance between two clusters according to the `Median link` merging rule.
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = \frac{d(C_i,C_k)}{2} +  \frac{d(C_j,C_k)}{2} - \frac{d(C_i,C_j)}{4}
/// ```
pub fn median_link(_size_i: usize, _size_j: usize, _size_k: usize,
                    distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {

    return 0.5 * distance_ik + 0.5 * distance_jk - 0.25 * distance_ij;
}

/// Calculates the distance between two clusters according to the `Centroid link` merging rule.
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = \frac{|C_i|}{|C_i|+|C_j|}d(C_i,C_k) +  \frac{|C_j|}{|C_i|+|C_j|}d(C_j,C_k) - \frac{|C_j||C_j|}{(|C_i|+|C_j|)^2}d(C_j,C_k)
/// ```
pub fn centroid_link(size_i: usize, size_j: usize, _size_k: usize,
                   distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {

    let d = 1.0 / ((size_i + size_j) as f32);
    return d * size_i as f32 * distance_ik + d * size_j as f32 * distance_jk - size_i as f32 * size_j as f32 * d * d * distance_ij;
}

/// Calculates the distance between two clusters according to the `minimum variance` merging rule (Ward's method).
///
/// The distance between any cluster $C_k$ and a new cluster created by merging $C_i$ and $C_j$ is defined as:
/// ```math
/// d(C_{i+j},C_k) = \frac{|C_i|+|C_k|}{|C_i|+|C_j|+|C_k|}d(C_i,C_k) +  \frac{|C_j|+|C_k|}{|C_i|+|C_j|+|C_k|}d(C_j,C_k) - \frac{|C_k|}{(|C_i|+|C_j|+|C_k|)}d(C_j,C_k)
/// ```
pub fn wards_method(size_i: usize, size_j: usize, size_k: usize,
                     distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {

    let d = 1.0 / ((size_i + size_j + size_k) as f32);
    return d * (size_i + size_k) as f32 * distance_ik + d * (size_j + size_k) as f32 * distance_jk - size_k as f32 * d * distance_ij;
}