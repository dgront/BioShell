//! Calculate sequence alignment. 
//!
//! This module provides implementation of Needleman-Wunsh and Smith-Waterman algorithms
//! for optimal sequence alignment with affine gap penalty
//!
//! # Global sequence alignment
mod global;
mod alignment_path;
mod alignment_reporter;
mod alignment_statistics;

pub use global::*;
pub use alignment_path::*;
pub use alignment_reporter::*;
pub use alignment_statistics::*;

