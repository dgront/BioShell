// functions and structs used by more than one clustering approach
mod points_set;
mod distance;

// re-export symbols to the top-most level of the module's name space
pub use points_set::{PointsWithDistance, EuclideanPoints};
pub use distance::{euclidean_distance_squared};

// each clustering method is placed in its own module
pub mod kd_tree;
pub mod optics;
pub mod em;




