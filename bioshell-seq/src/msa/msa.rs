use std::io::BufRead;
use std::collections::{HashMap, HashSet};

use crate::errors::SequenceError;
use crate::sequence::{FastaIterator, Sequence, StockholmIterator};

#[derive(Default, Clone, Debug)]
/// Multiple Sequence Alignment is a `Vec` of  [`Sequence`](Sequence) objects of the same length.
///
pub struct MSA {
    /// a sequence is represented as a vector of u8 bytes
    msa: Vec<Sequence>,
}

impl MSA {

    pub fn from_stockholm_reader<R: BufRead>(reader: &mut R) -> Result<Self, SequenceError> {
        let msa = StockholmIterator::new( reader).into_iter().collect();
        Self::check_msa(&msa)?;
        return Ok(MSA { msa });
    }

    pub fn from_fasta_reader<R: BufRead>(reader: &mut R) -> Result<Self, SequenceError> {
        let msa = FastaIterator::new( reader).into_iter().collect();
        Self::check_msa(&msa)?;
        return Ok(MSA { msa });
    }

    pub fn len(&self) ->usize { self.msa[0].len() }

    pub fn n_seq(&self) ->usize { self.msa.len() }

    fn check_msa(msa: &Vec<Sequence>) -> Result<(), SequenceError> {
        // get unique sequence length values
        let unique_len: HashSet<usize> = msa.iter().map(|s| s.len()).collect();
        // It's OK when the set contains just a single value
        if unique_len.len()==1 { return Ok(()) }
        // Otherwise, find the most popular length and report an error
        let mut count_by_len: HashMap<usize, usize> = unique_len.iter()
            .map(|&c| (c, 0))
            .collect();
        for s in msa {
            *count_by_len.get_mut(&s.len()).unwrap() += 1;
        }
        let mut count_vec: Vec<(&usize, &usize)> = count_by_len.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));

        return Err(SequenceError::AlignedSequencesOfDifferentLengths {
            length_expected: *count_vec[0].0,
            length_found: *count_vec[1].0,
        });
    }
}