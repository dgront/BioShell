//! Hierarchical agglomerate clustering algorithm
//!
mod hierarchical;
mod distance_matrix;
mod merging_rules;

pub use hierarchical::*;
pub use merging_rules::*;
use distance_matrix::*;

