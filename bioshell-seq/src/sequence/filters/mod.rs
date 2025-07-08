//! Sequence filters are be used to select sequences from a larger pool of sequences.
//!
//! Each [`SequenceFilter`] returns `true` or `false` for a given sequence object.
//! Sequences may be filtered by length, description, amino acid composition, sequence motif, and other criteria.
//! For example, [`ContainsAA`] filter returns `true` if a given [`Sequence`] contains a given amino acid type at least `n` times
//!
//! # Examples
//! ```rust
//! use bioshell_seq::sequence::Sequence;
//! use bioshell_seq::sequence::filters::{SequenceFilter, ContainsAA};
//! let fdx_seq = Sequence::from_str("test_seq", "VVFGCKRCGKCRDVCPVGAIYEENELAKIDTEKCNLCMKCIDECTNRSIIYME");
//! let filter = ContainsAA{res_type: 'C',min_cnt: 7};
//! assert_eq!(filter.filter(&fdx_seq),true);
//! ```

mod sequence_filters;
pub use sequence_filters::*;
use crate::sequence::Sequence;