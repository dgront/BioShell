//! Efficient and clean library for processing biomacromolecular structures.
//!

#![allow(clippy::needless_return)]
pub mod nerf;

mod internal_coordinates_definitions;
mod kinematic_atom_tree;
mod builder_errors;

pub use internal_coordinates_definitions::{RelaviveResidueLocator, InternalAtomDefinition, InternalCoordinatesDatabase};
pub use kinematic_atom_tree::{KinematicAtomTree};
pub use builder_errors::*;