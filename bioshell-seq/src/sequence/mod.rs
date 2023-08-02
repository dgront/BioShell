//! Provides [`Sequence`](crate::sequence::Sequence) and [`SequenceProfile`](crate::sequence::SequenceProfile) stucts
//! as well as mumerous utility functions to operate on them


mod sequence;
mod sequence_profile;
mod residue_type_mapping;
mod sequence_filters;

pub use sequence::*;
pub use sequence_filters::*;
pub use sequence_profile::*;
pub use residue_type_mapping::*;