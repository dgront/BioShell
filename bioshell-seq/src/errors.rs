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
}

#[allow(unused_imports)]
#[allow(dead_code)]
mod test_error {
    #[allow(unused_imports)]
    use std::io::BufReader;
    use crate::errors::SequenceError;
    use crate::msa::MSA;

    #[allow(non_upper_case_globals)]
    static fasta: &'static str = "> seq A
MTYKL
> seq B
MTYK-
> seq C
MTYK";

    #[test]
    fn read_msa_fasta() {
        let mut fas_reader = BufReader::new(fasta.as_bytes());
        let msa_error = MSA::from_fasta_reader(&mut fas_reader).err().unwrap();
        let expected_err = SequenceError::AlignedSequencesOfDifferentLengths {
            length_expected: 5,
            length_found: 4,
        };
        assert_eq!(msa_error, expected_err);
    }
}
