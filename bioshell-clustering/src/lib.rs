// functions and structs used by more than one clustering approach
mod points_set;
mod distance;

// re-export symbols to the top-most level of the module's name space
pub use points_set::{DistanceByIndex, CartesianPoints, Neighbor, NeighborsOf};
pub use distance::euclidean_distance_squared;

// each clustering method is placed in its own module
pub mod optics;
pub mod em;
pub mod kmeans;




