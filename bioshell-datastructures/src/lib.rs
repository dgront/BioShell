// since the crate is quite small so far, all its structs are in the root namespace
mod tree;
// re-export symbols to the top-most level of the module's name space
pub use tree::*;

// k-d tree are in a separate crate to keep it easy to follow
pub mod kd_tree;

mod distance;
pub use distance::euclidean_distance_squared;





