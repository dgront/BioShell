mod optics;
mod em;
pub mod kd_tree;
mod points_set;

pub use points_set::{PointsWithDistance, EuclideanPoints};

pub use optics::{Optics, Neighbor, NeighborsOf, OpticsPoints};
pub use em::{expectation_maximization};


