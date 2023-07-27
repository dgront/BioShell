use core::fmt::{Formatter, Display};
use crate::msa::MSA;

use crate::sequence::{ProfileColumnOrder};

#[derive(Clone, Debug)]
/// Represents a sequence profile.
///
/// Sequence profile describes preference for each of the 20 standard amino acid residue types
/// (possibly also with the gap symbol) at each of the residue positions in the query sequence.
/// For a protein sequence of N residues, these preferences are stored as a `Nx20` probability matrix
/// normalised row-wise, i.e. they sum up to `1.0` for each sequence position. These probabilities
/// reflect how likely is to find any amino acid type (say, GLY or TYR) at that position.
///
/// Integer indexes are used to refer to both a residue position in a sequence and an amino acid (or nucleotide) type;
/// both indexes start from 0. Residue type indexes must be consistent with the [`ProfileColumnOrder`]
/// object used by a given profile.
pub struct SequenceProfile {
    mapping: ProfileColumnOrder,
    data: Vec<Vec<f32>>,
}

impl SequenceProfile {
    /// Creates a sequence profile from a given  [`MSA`] (Multiple Sequence Alignment) object
    pub fn new(mapping: ProfileColumnOrder, msa: &MSA) -> SequenceProfile {

        let n = msa.len();
        let k = mapping.size();
        let mut data: Vec<Vec<f32>> = vec![vec![0.0; k]; n];

        for sequence in msa.sequences() {
            for idx in 0..n {
                let aa_idx = mapping.type_to_index(&sequence.u8(idx)) as usize;
                data[idx][aa_idx] += 1.0;
            }
        }

        for idx in 0..n {
            let sum: f32 = data[idx].iter().sum();
            for aa_idx in 0..k {
                data[idx][aa_idx] /= sum;
            }
        }

        SequenceProfile{mapping, data}
    }

    /// Returns a probability for a given sequence position and a residue type.
    pub fn fraction(&self, pos: usize, aa: usize) -> f32 { self.data[pos][aa] }

    /// Provides immutable access to the residue ordering used by this sequence profile
    pub fn column_order(&self) -> &ProfileColumnOrder { &self.mapping }

    /// Says which column holds statistics for an amino acids given by its single-char code
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::{ProfileColumnOrder, Sequence, SequenceProfile};
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "actg")]).unwrap();
    /// let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(), &msa);
    ///
    /// assert_eq!(profile.column_index(&'a'), 0);
    /// assert_eq!(profile.column_index(&'g'), 2);
    /// ```
    pub fn column_index(&self, aa: &char) -> usize { self.mapping.letter_to_index(aa) }

    /// Returns the number of sequence positions in this sequence profile
    pub fn len(&self) -> usize { self.data.len() }
}

impl Display for SequenceProfile {
    /// Prints a sequence profile as a nice table
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (i, row) in self.data.iter().enumerate() {
            write!(f, "{i:4} ").ok();
            for val in row {
                write!(f, "{val:5.3} ").ok();
            }
            writeln!(f, "").ok();
        }
        Ok(())
    }
}