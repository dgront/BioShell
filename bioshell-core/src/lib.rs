//! Core structs and utility functions used by bioshell crates
//!

#![allow(clippy::needless_return)]

pub mod io;
mod matrix;
pub use matrix::*;

mod vec3;
pub use vec3::*;

mod has_cartesians;
mod assertions;

pub use has_cartesians::HasCartesians;