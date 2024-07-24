use core::fmt::{Formatter, Display};
use std::string::String;

use crate::msa::MSA;

use crate::sequence::{ProfileColumnOrder, Sequence};

#[derive(Clone, Debug)]
/// Represents a sequence profile.
///
/// A protein sequence profile describes preference for each of the 20 standard amino acid residue types
/// (possibly also with the gap symbol) at each of the residue positions in the query sequence.
/// For a protein sequence of N residues, these preferences are stored as a `Nx20` probability matrix
/// normalised row-wise, i.e. they sum up to `1.0` for each sequence position. These probabilities
/// reflect how likely is to find any amino acid type (say, GLY or TYR) at that position in a given protein family.
/// Similiarily, such a sequence profile may be constructed from DNA or RNA sequences.
///
/// Integer indexes are used to refer to both a residue position in a sequence and an amino acid (or nucleotide) type;
/// both indexes start from 0. Residue type indexes must be consistent with the [`ProfileColumnOrder`]
/// object used by a given profile.
///
/// # Example
/// ```rust
/// # use bioshell_seq::sequence::{Sequence, SequenceProfile, ProfileColumnOrder};
/// # use bioshell_seq::msa::MSA;
/// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "cttaga"),
///                                    Sequence::from_str("seq-2", "ctcaga"),
///                                    Sequence::from_str("seq-3", "ctaagg"),
///                                    Sequence::from_str("seq-4", "cttaga")]).unwrap();
/// let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(), &msa);
/// assert_eq!(profile.most_probable_string(), String::from("cttaga"));
/// ```
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
    pub fn fraction(&self, pos: usize, column: usize) -> f32 { self.data[pos][column] }

    /// Provides immutable access to the residue ordering used by this sequence profile
    pub fn column_order(&self) -> &ProfileColumnOrder { &self.mapping }

    /// Says which column holds statistics for a residue type given by its single-char code
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::{ProfileColumnOrder, Sequence, SequenceProfile};
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "actg")]).unwrap();
    /// // Column order of this profile is: 'a', 'c', 'g' and 't', as defined by ProfileColumnOrder::dna_standard()
    /// let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(), &msa);
    ///
    /// assert_eq!(profile.column_index(&'a'), 0);  // 'a' has index 0
    /// assert_eq!(profile.column_index(&'g'), 2);  // 'g' has index 2
    /// ```
    pub fn column_index(&self, aa: &char) -> usize { self.mapping.letter_to_index(aa) }

    /// Character that represents the residue type corresponding to a given profile column
    ///
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::{ProfileColumnOrder, Sequence, SequenceProfile};
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "actg")]).unwrap();
    /// // Column order of this profile is: 'a', 'c', 'g' and 't', as defined by ProfileColumnOrder::dna_standard()
    /// let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(), &msa);
    ///
    /// assert_eq!(profile.column_char(0u8), 'a');  // 'a' has index 0
    /// assert_eq!(profile.column_char(2u8), 'g');  // 'g' has index 2
    /// ```
    pub fn column_char(&self, column: u8) -> char { self.mapping.index_to_letter(column) }

    /// Returns the number of sequence positions in this sequence profile
    pub fn len(&self) -> usize { self.data.len() }

    /// Returns the index of a column that holds maximum probability for a given profile position
    pub fn most_probable_res(&self, pos: usize) -> usize {
        self.data[pos].iter().enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(index, _)| index).unwrap()
    }

    /// Returns the most probable sequence.
    ///
    /// Returned string contains characters that are the occur most often at each position of this profile
    /// # Example
    /// ```rust
    /// # use bioshell_seq::sequence::{Sequence, SequenceProfile, ProfileColumnOrder};
    /// # use bioshell_seq::msa::MSA;
    /// let msa = MSA::from_sequences(vec![Sequence::from_str("seq-1", "cttaga"),
    ///                                    Sequence::from_str("seq-2", "ctcaga"),
    ///                                    Sequence::from_str("seq-3", "ctaagg"),
    ///                                    Sequence::from_str("seq-4", "cttaga")]).unwrap();
    /// let profile = SequenceProfile::new(ProfileColumnOrder::dna_standard(), &msa);
    /// assert_eq!(profile.most_probable_sequence().to_string(0), String::from("cttaga"));
    /// ```
    pub fn most_probable_sequence(&self) -> Sequence {
        return Sequence::new(&String::from("most probable sequence"), &self.most_probable_string());
    }

    /// Returns the most probable sequence as a string.
    ///
    /// This methods works as [`most_probable_sequence()`](Self::most_probable_sequence()), but returns just a string representation
    /// of a sequence. Use this method when a full Sequence is not needed.
    pub fn most_probable_string(&self) -> String {
        let mut resids = vec![];
        for i in 0..self.len() {
            resids.push(self.column_char(self.most_probable_res(i) as u8));
        }
        let s: String = resids.into_iter().collect();
        return s;
    }
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