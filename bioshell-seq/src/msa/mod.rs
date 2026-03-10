//! Provides [`MSA`] structs and functions that operate on such alignments
mod msa;
pub use msa::*;

mod display_msa;

mod stockholm_io;
mod stockholm_msa;
pub use stockholm_msa::*;

