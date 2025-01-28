//! Hierarchical agglomerate clustering algorithm
//!
//! In this algorithm, the clustering procedure builds a binary tree, or dendrogram, of clusters.
//! At each step, the two clusters that are closest to each other are merged into a new cluster.
//! The algorithm stops when all clusters have been merged into one big cluster containing all objects.
//! An example of a dendrogram resulting from clustering the following real values:
//! ``data = vec![6.0, 2.0, 2.5, 5.9]`` is shown below:
//!
//! ```text
//!         6
//!       /   \
//!      4     5
//!     / \   / \
//!    0   3  1  2
//! ```
//!
//! where numbers ``0, 3, 1, 2`` in the bottom row represent the clustering data items, e.g. ``data[0] = 6.0``.
//!
//! This implementation does not use the data points explicitly. Instead, it uses
//! a distance function to compute the distance between two object identified by their index.
//! For example, to cluster the following real values: ``[6.0, 2.0, 2.5, 5.9]``, we define the distance function as:
//! ```
//! let data: Vec<f32> = vec![6.0, 2.0, 2.5, 5.9];
//! let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();
//! # assert_eq!(distance_fn(0, 1), 4.0);
//! ```
//!
//! Now we can cluster the data using the hierarchical clustering algorithm:
//! ```
//! use bioshell_clustering::hierarchical::hierarchical_clustering;
//! use bioshell_clustering::hierarchical::strategies::single_link;
//! # let data: Vec<f32> = vec![6.0, 2.0, 2.5, 5.9];
//! # let distance_fn = |ia: usize, ib: usize| (data[ia] - data[ib]).abs();
//! let mut clustering = hierarchical_clustering(data.len(), distance_fn, &single_link);
//! ```
//! which returns the root node of a binary clustering tree.
//!
mod hierarchical;
mod clustering_matrix;
pub mod strategies;
mod distance_matrix;

pub use hierarchical::*;
pub use distance_matrix::DistanceMatrix;
use clustering_matrix::*;

