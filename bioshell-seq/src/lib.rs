//! Efficiently process biological sequences and their annotations.
//!
//! A brief summary of the `bioshell-seq` library features are provided below. Documentation
//! of command line executables may be found in [cookbook](documentation)
//!
//! This crate is designed to facilitate common tasks related to protein and nucleic sequences
//! as well as sequence profiles. It provides:
//!   - [`Sequence`](sequence::Sequence) struct to store an amino acid or a nucleic sequence,
//!     which can be loaded e.g. from the [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
//!     format with [`FastaIterator`](sequence::FastaIterator)
//!   - [`MSA`](msa::MSA) struct that stores a number of sequences as a `Vec<Sequence>`,
//!     which are assumed to be aligned with one another. [`MSA`](msa::MSA) struct data
//!     can be loaded from the following file formats:
//!     - [FASTA](https://en.wikipedia.org/wiki/FASTA_format), including the ``a3m`` variant
//!     - [Stockholm](https://sonnhammer.sbc.su.se/Stockholm.html)
//!   - [`SequenceRecord`](crate::sequence::SequenceRecord) struct to store the full entry of a sequence and its annotations.
//!     It can be loaded from a:
//!     - [GeneBank](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) file
//!     - [UniProt](https://web.expasy.org/docs/userman.html) file
//!   - [`SequenceProfile`](crate::sequence::SequenceProfile) struct to store a sequence profile
//!     with a given [`ProfileColumnOrder`](crate::sequence::ProfileColumnOrder). A [`SequenceProfile`](crate::sequence::SequenceProfile)
//!     object may be created from an [`MSA`](crate::msa::MSA) instance
//!
//!
pub mod documentation;

mod errors;
pub mod chemical;
pub mod sequence;
pub mod scoring;
pub mod alignment;
pub mod msa;

pub use errors::{SequenceError, ScoringError};
