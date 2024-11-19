//! Provides merging rules for hierarchical clustering.
//!
//! The merging rules are used to determine the distance between two clusters, based on the
//! distances between the clusters elements themselves. There are four merging rules implemented:
//!
//!  - `Single link`: the distance between two clusters is the **minimum** distance between any two elements of the two clusters.
//!  - `Complete link`: the distance between two clusters is the **maximum** distance between any two elements of the two clusters.
//!  - `Average link`: the distance between two clusters is the **average** distance between any two elements of the two clusters.
//!  - `Centroid link`: the distance between two clusters is the distance between the centroids of the two clusters.
//!
//! These merging rules can be evaluated in constant time, which requires some additional information.
//! Specifically, the distance between any cluster `L` and a cluster `K`  just created by merging
//! two clusters `I` and `J`,  requires the sizes of the clusters `I`, `J` and `K`,
//! and the distances between the clusters `I` and `J`, `I` and `K`, and `J` and `K`.
//!

/// Calculates the distance between two clusters according to the Single link merging rule.
pub fn single_link(size_i: usize, size_j: usize, size_k: usize,
                   distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {
    distance_ik.min(distance_jk)
}

/// Calculates the distance between two clusters according to the `Complete link` merging rule.
pub fn complete_link(size_i: usize, size_j: usize, size_k: usize,
                   distance_ij: f32, distance_ik: f32, distance_jk: f32) -> f32 {
    distance_ik.max(distance_jk)
}