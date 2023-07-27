use std::io::BufRead;
use std::collections::{HashMap, HashSet};

use crate::errors::SequenceError;
use crate::sequence::{FastaIterator, Sequence, StockholmIterator};

#[derive(Default, Clone, Debug)]
/// Multiple Sequence Alignment is a `Vec` of  [`Sequence`](Sequence) objects of the same length.
///
pub struct MSA {
    /// holds all sequences aligned
    msa: Vec<Sequence>,
}

impl MSA {

    /// Create an MSA by consuming the given sequences
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::sequence::Sequence;
    /// use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
    ///                                    Sequence::from_str("seq-2", "PERV")]).unwrap();
    /// assert_eq!(msa.n_seq(), 2);
    /// ```
    pub fn from_sequences(sequences: Vec<Sequence>) -> Result<Self, SequenceError> {
        let msa = sequences;
        Self::check_msa(&msa)?;
        return Ok(MSA { msa });
    }

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

    /// Returns the length of each sequence of this alignment
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::Sequence;
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
    ///                                    Sequence::from_str("seq-2", "PERV")]).unwrap();
    /// assert_eq!(msa.len(), 4);
    /// ```
    pub fn len(&self) ->usize { self.msa[0].len() }

    /// Returns the number of sequences in this alignment
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::Sequence;
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
    ///                                    Sequence::from_str("seq-2", "PERV")]).unwrap();
    /// assert_eq!(msa.n_seq(), 2);
    /// ```
    pub fn n_seq(&self) ->usize { self.msa.len() }

    /// Provide immutable access to the sequences of this alignment
    pub fn sequences(&self) -> &Vec<Sequence> { &self.msa }

    /// Iterate over a given column of this MSA
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::Sequence;
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
    ///                                    Sequence::from_str("seq-2", "P-RV"),
    ///                                    Sequence::from_str("seq-3", "PDRV")]).unwrap();
    /// let column1: Vec<u8> = msa.nth_column_iter(1).collect();
    /// assert_eq!(column1, [b'E', b'-', b'D']);
    /// ```
    pub fn nth_column_iter(&self, n: usize) -> IteratorOverColumn {

        IteratorOverColumn{
            msa: &self.msa,
            index: 0,
            column: n
        }
    }

    /// Counts identical residues between two sequences of this MSA.
    ///
    /// Matching gap symbols are not included in the count.
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::Sequence;
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "PERF"),
    ///                                    Sequence::from_str("seq-2", "P-RV"),
    ///                                    Sequence::from_str("seq-3", "P_RV")]).unwrap();
    /// assert_eq!(msa.identical_count(1,2), 3);
    /// assert_eq!(msa.identical_count(0,1), 2);
    /// ```
    pub fn identical_count(&self, i_seq: usize, j_seq: usize) -> usize {
        let mut ret = 0;
        let si = &self.msa[i_seq];
        let sj = &self.msa[j_seq];
        for i in 0..self.len() {
            if si.u8(i)==sj.u8(i) && si.u8(i)!=b'-' && si.u8(i)!=b'_' {
                ret += 1;
            }
        }
        return ret;
    }

    /// Computes the sequence identity fraction between two sequences of this MSA.
    ///
    /// This methods merely returns the `identical_count()` divided by the length of the sequences
    pub fn identical_fraction(&self, i_seq: usize, j_seq: usize) -> f64 {

        return self.identical_count(i_seq, j_seq) as f64 / self.len() as f64;
    }

    fn check_msa(msa: &Vec<Sequence>) -> Result<(), SequenceError> {
        // get unique sequence descriptions
        let unique_desc: HashSet<&str> = msa.iter().map(|s| s.description()).collect();
        // If we have fewer description than sequences - there must be duplicates
        if unique_desc.len() < msa.len() {
            let mut count_by_desc: HashMap<&str, usize> = unique_desc.iter()
                .map(|&c| (c, 0)).collect();
            for s in msa {
                *count_by_desc.get_mut(&s.description()).unwrap() += 1;
            }
            let mut count_vec: Vec<(&&str, &usize)> = count_by_desc.iter().collect();
            count_vec.sort_by(|a, b| b.1.cmp(a.1));
            return Err(SequenceError::IdenticalSequenceDescriptions {
                description: String::from(*count_vec[0].0),
            });
        }

        // get unique sequence length values
        let unique_len: HashSet<usize> = msa.iter().map(|s| s.len()).collect();
        // If the set contains more than just a single value, find the most popular length and report an error
        if unique_len.len() != 1 {
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
        return Ok(());
    }
}

pub struct IteratorOverColumn<'a> {
    msa: &'a Vec<Sequence>,
    index: usize,
    column: usize,
}

impl Iterator for IteratorOverColumn<'_> {
    type Item = u8;
    fn next(&mut self) -> Option<u8> {
        let result = if self.index >= self.msa.len() { None }
                        else { Some(*&self.msa[*(&self.index)].u8(*(&self.column))) };
        self.index += 1;
        return result;
    }
}