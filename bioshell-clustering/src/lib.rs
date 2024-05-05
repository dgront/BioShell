// functions and structs used by more than one clustering approach
mod points_set;

// re-export symbols to the top-most level of the module's name space
pub use points_set::{DistanceByIndex, CartesianPoints, Neighbor, NeighborsOf};

// each clustering method is placed in its own module
pub mod optics;
pub mod em;
pub mod kmeans;




