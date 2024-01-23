//! Builds 3D biomacromolecular structures.
//!
//! ``bioshell-builder`` crate builds a bimacromolecular structure (a polypeptide or a nucleic acid chain)
//! from its sequence.
//!
//! ```
//! use bioshell_builder::InternalCoordinatesDatabase;
//! let library = InternalCoordinatesDatabase::from_cif_directory("./data/")
//!     .expect("Can't load residue definitions from a given folder: './data/'");
//! # assert!(library.count_definitions() > 0);
//!
//! ```


#![allow(clippy::needless_return)]
pub mod nerf;

mod internal_coordinates_definitions;
mod kinematic_atom_tree;
mod builder_errors;

pub use internal_coordinates_definitions::{RelaviveResidueLocator, InternalAtomDefinition, InternalCoordinatesDatabase};
pub use kinematic_atom_tree::{KinematicAtomTree};
pub use builder_errors::*;