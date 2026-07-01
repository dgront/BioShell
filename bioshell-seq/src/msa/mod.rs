//! Provides [`MSA`] structs and functions that operate on such alignments
mod msa;
pub use msa::*;

mod display_msa;

mod stockholm_io;
// pub use stockholm_io::*;

mod stockholm_msa;
pub use stockholm_msa::*;

mod msa_calculations;
mod clustalw_io;

pub use msa_calculations::*;
