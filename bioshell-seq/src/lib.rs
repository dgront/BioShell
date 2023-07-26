//! Efficiently process biological sequences.
//!
//! This crate is designed to facilitate common tasks related to protein and nucleic sequences
//! as well as sequence profiles. It provides:
//!  - [`Sequence`](Sequence) struct to store an amino acid or a nucleic sequence, which can be loaded from the following file formats:
//!    - [FASTA](https://en.wikipedia.org/wiki/FASTA_format), including the ``a3m`` variant
//!    - [Stockholm](https://sonnhammer.sbc.su.se/Stockholm.html)
//!

pub mod chemical;
pub mod sequence;
pub mod msa;
