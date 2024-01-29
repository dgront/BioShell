//! Calculate sequence alignment scores
//!
//! # Substitution matrices
//! [SubstitutionMatrix] struct holds substitution (similarity) matrix scores; it provides a score
//! for aligning any two amino acid types. The [`bioshell-seq`](crate) crate defines a few most popular matrices,
//! listed in the [SubstitutionMatrixList] enum. Such a matrix may be loaded as shown below:
//!
//! ```
//! use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
//! let blosum80 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM80);
//! ```
//!
//! This data can be also loaded from data in the NCBI format, e.g:
//!
//! ```
//! const BLOSUM62: &str = "# Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0.
//!    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
//! A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -1 -1 -4
//! R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1 -2  0 -1 -4
//! N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  4 -3  0 -1 -4
//! D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4 -3  1 -1 -4
//! C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -1 -3 -1 -4
//! Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0 -2  4 -1 -4
//! E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1 -3  4 -1 -4
//! G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -4 -2 -1 -4
//! H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0 -3  0 -1 -4
//! I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3  3 -3 -1 -4
//! L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4  3 -3 -1 -4
//! K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0 -3  1 -1 -4
//! M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3  2 -1 -1 -4
//! F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3  0 -3 -1 -4
//! P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -3 -1 -1 -4
//! S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0 -2  0 -1 -4
//! T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1 -1 -1 -4
//! W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -2 -2 -1 -4
//! Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -1 -2 -1 -4
//! V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3  2 -2 -1 -4
//! B -2 -1  4  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4 -3  0 -1 -4
//! J -1 -2 -3 -3 -1 -2 -3 -4 -3  3  3 -3  2  0 -3 -2 -1 -2 -1  2 -3  3 -3 -1 -4
//! Z -1  0  0  1 -3  4  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -2 -2 -2  0 -3  4 -1 -4
//! X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4
//! * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1";
//!
//! use std::io::BufReader;
//! use bioshell_seq::scoring::SubstitutionMatrix;
//! let blosum62 = SubstitutionMatrix::ncbi_matrix_from_buffer(BufReader::new(BLOSUM62.as_bytes())).unwrap();
//! assert_eq!(blosum62.score_by_aa('C' as u8, 'C' as u8), 9);
//! ```
//!
//! # Scoring
//! A substitution matrix may be used to score similarity between two amino acid types. The
//! [SubstitutionMatrix] struct API allows access the score either by amino acid one-character codes
//! (an ``u8`` byte representing a respective character) or by their internal indexes (``u8`` values
//! from 0 to 20, both inclusive, with the index 20 represents the ``'X'`` amino acid)
//!
//! ```
//! use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
//! let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
//! // --- Amino acid similarity scores can be accessed directly by one-letter codes
//! assert_eq!(blosum62.score_by_aa('W' as u8, 'W' as u8), 11);
//! assert_eq!(blosum62.score_by_aa('X' as u8, 'X' as u8), -1);
//! // --- Alternatively, one ca cache internal amino acid indexes for direct access
//! let ala_idx = blosum62.aa_index('A' as u8);
//! let trp_idx = blosum62.aa_index('W' as u8);
//! assert_eq!(blosum62.score_by_index(ala_idx, trp_idx), -3);
//! // let's double check if the matrix is symmetric
//! assert_eq!(blosum62.score_by_index(trp_idx, ala_idx), -3);
//! ```
mod substitution_matrix;
mod similarity_score;

pub use substitution_matrix::*;
pub use similarity_score::*;

