use displaydoc::Display;
use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
/// Errors that may be thrown while processing sequences or multiple sequence alignments
pub enum SequenceError {
    #[error("Expected aligned sequences of length {length_expected}, found sequences of lengths: {length_found}")]
    /// A sequence has been found that doesn't match by length to the alignment
    AlignedSequencesOfDifferentLengths {
        /// Expected length
        length_expected: usize,
        /// Found length
        length_found: usize,
    },
    #[error("Expected the new aligned sequence of length {length_expected}, found sequence of length {length_found}")]
    /// Attempted to add a sequence that doesn't match by length to the alignment
    NewSequenceOfDifferentLength {
        /// Expected length
        length_expected: usize,
        /// Found length
        length_found: usize,
    },

    #[error("The following description: '{description}' has been found in more than one sequence")]
    /// The given `description` has been found in more than one sequence
    IdenticalSequenceDescriptions {
        /// multiplied description
        description: String,
    },

    #[error("Can't locate a sequence for the given ID: '{seq_id}'")]
    /// Can't locate a sequence for the given `seq_id`
    InvalidSequenceID {
        /// the ID of the missing sequence
        seq_id: String,
    },

    #[error("Invalid Stockholm file format")]
    /// Invalid Stockholm file format
    InvalidStockholmFormat {
        /// the line causing the error
        line: String,
        /// error description
        description: String,
    },

    #[error("Invalid Fasta file format")]
    /// Invalid Stockholm file format
    InvalidFastaFormat {
        /// the line causing the error
        line: String,
        /// error description
        description: String,
    },

    #[error("General I/O error occurred while reading or writing a sequence file")]
    /// I/O error occurred while reading a sequence file
    Io(#[from] std::io::Error),

    #[error("The provided string is not a valid sequence pattern: {pattern}")]
    /// The provided string is not a valid sequence pattern
    IncorrectSequencePattern{ pattern: String },


    #[error("The provided set of sequences or a multiple sequence alignment is empty")]
    /// The provided set of sequences or a multiple sequence alignment is empty
    NoInputSequences,
}

#[derive(Debug, Error, Display, PartialEq)]
#[non_exhaustive]
/// Errors that may be thrown while loading or using a substitution matrix
pub enum ScoringError {
    /// The file: {file_name} can't be opened for reading
    FileNotFound {
        /// name of the missing file
        file_name: String,
    },
    /// Reading occurred while reading a substitution matrix
    ReadingError,
    /// Given amino acid index: {given_index} is too high, maximum value is 20
    AminoAcidIndexTooHigh {
        /// Expected length
        given_index: u8,
    },
    /// The following line of a NCBI matrix file is not formatted correctly: {line}
    IncorrectNCBIFormat {
        /// the incorrectly formatted line that broke the code
        line: String,
    },
    /// The following entry: {value} found in line can't be parsed to i32 type; the problematic line was: {line}
    CantParseNCBIEntry {
        /// the incorrectly formatted line that broke the code
        line: String,
        value: String
    },
}

#[allow(unused_imports)]
#[allow(dead_code)]
mod test_error {
    #[allow(unused_imports)]
    use std::io::BufReader;
    use crate::errors::SequenceError;
    use crate::msa::MSA;
    use crate::sequence::Sequence;
    use crate::SequenceError::{IdenticalSequenceDescriptions, AlignedSequencesOfDifferentLengths};

    #[allow(unused_macros)]
    macro_rules! assert_err {
        ($expression:expr, $($pattern:tt)+) => {
          match $expression {
                $($pattern)+ => (),
                ref e => panic!("expected `{}` but got `{:?}`", stringify!($($pattern)+), e),
            }
        }
    }
    #[test]
    fn check_different_lengths() {
        let msa_error = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                                 Sequence::from_str("seq-2", "P-RV"),
                                                 Sequence::from_str("seq-3", "PRV")]);
        assert_err!(msa_error, Err(AlignedSequencesOfDifferentLengths{length_expected: e, length_found: f}) if e==4 && f == 3);
    }

    #[test]
    fn check_identical_descriptions() {
        let msa_error = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                        Sequence::from_str("seq-1", "P-RV")]);
        assert_err!(msa_error, Err(IdenticalSequenceDescriptions{description: desc}) if desc == "seq-1".to_string());
    }
}
