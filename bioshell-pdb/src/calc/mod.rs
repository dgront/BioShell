//! Functions that calculate various geometric properties such as distances, angles etc.
mod simple_geometric;
mod vec3;

pub use simple_geometric::*;
pub use vec3::{Vec3, planar_angle2, planar_angle3, dihedral_angle4,
               random_unit_versor, random_point_nearby};