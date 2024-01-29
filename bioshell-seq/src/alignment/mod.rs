//! Calculate sequence alignment. 
//!
//! This module provides implementation of Needleman-Wunsh and Smith-Waterman algorithms
//! for optimal sequence alignment with affine gap penalty
//!
//! # Global sequence alignment
mod global;
mod aligner;

pub use global::*;
pub use aligner::*;

