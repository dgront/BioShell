use displaydoc::Display;
use thiserror::Error;

#[derive(Debug, Error, Display, PartialEq)]
#[non_exhaustive]
/// Errors that may be thrown by `bioshell-seq` crate
pub enum SequenceError {
    /// Expected aligned sequences of length {length_expected}, found sequences of lengths: {length_found}
    AlignedSequencesOfDifferentLengths {
        /// Expected length
        length_expected: usize,
        /// Found length
        length_found: usize,
    },
    /// Expected new aligned sequence of length {length_expected}, found sequence of length {length_found}
    NewSequenceOfDifferentLength {
        /// Expected length
        length_expected: usize,
        /// Found length
        length_found: usize,
    },
    /// The following description: "{description}" has been found in more than one sequence
    IdenticalSequenceDescriptions {
        /// multiplied description
        description: String,
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

    #[test]
    fn check_different_lengths() {
        let msa_error = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                                 Sequence::from_str("seq-2", "P-RV"),
                                                 Sequence::from_str("seq-3", "PRV")]).err().unwrap();
        let expected_err = SequenceError::AlignedSequencesOfDifferentLengths {
            length_expected: 4,
            length_found: 3,
        };
        assert_eq!(msa_error, expected_err);
    }

    #[test]
    fn check_identical_descriptions() {
        let msa_error = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
                                        Sequence::from_str("seq-1", "P-RV")]).err().unwrap();
        let expected_err = SequenceError::IdenticalSequenceDescriptions {
            description: "seq-1".to_string(),
        };
        assert_eq!(msa_error, expected_err);
    }
}
