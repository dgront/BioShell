use core::fmt::{Display, Formatter};

use crate::sequence::{ResidueTypeOrder, Sequence};

#[derive(Clone, Debug)]
/// Represents a sequence profile.
///
/// Sequence profile describes preference for each of the 20 standard amino acid residue types
/// (possibly also with the gap symbol) at each of the residue positions in the query sequence.
/// For a protein sequence of N residues, these preferences are stored as a `Nx20` probability matrix
/// normalised row-wise, i.e. they sum up to `1.0` for each sequence position. These probabilities
/// reflect how likely is to find any amino acid type (say, GLY or TYR) at that position.
///
/// Integer indexes are used to refer to both a residue position in a sequence and an amino acid (ro nucleotide) type;
/// both indexes start from 0. Residue type indexes must be consistent with the [`ResidueTypeOrder`]
/// object used by a given profile.
pub struct SequenceProfile {
    mapping: ResidueTypeOrder,
    data: Vec<Vec<f32>>,
}

impl SequenceProfile {
    /// Creates a sequence profile from a given pool of aligned sequences
    pub fn new(mapping: ResidueTypeOrder, msa: &Vec<Sequence>) -> SequenceProfile {
        let n = msa[0].len();
        let k = mapping.size();
        let mut data: Vec<Vec<f32>> = vec![vec![0.0; k]; n];

        for sequence in msa {
            if sequence.len() != n {
                eprintln!("\nSequence of incorrect length! Is: {}, should be: {}. The sequence skipped::\n {}\n",
                          sequence.len(), n, sequence);
                continue;
            }
            for idx in 0..n {
                let aa_idx = mapping.type_to_index(&sequence.seq()[idx]) as usize;
                data[idx][aa_idx] += 1.0;
            }
        }

        for idx in 0..n {
            let sum: f32 = data[idx].iter().sum();
            for aa_idx in 0..k {
                data[idx][aa_idx] /= sum;
            }
        }

        SequenceProfile { mapping, data }
    }

    /// Returns a probability for a given sequence position and a residue type.
    pub fn fraction(&self, pos: usize, aa: usize) -> f32 {
        self.data[pos][aa]
    }

    /// Provides immutable access to the residue ordering used by this sequence profile
    pub fn column_order(&self) -> &ResidueTypeOrder {
        &self.mapping
    }

    /// Returns the number of sequence positions in this sequence profile
    pub fn len(&self) -> usize {
        self.data.len()
    }
}

impl Display for SequenceProfile {
    /// Prints a sequence profile as a nice table
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (_i, row) in self.data.iter().enumerate() {
            //write!(f, "{i:4} ");
            for val in row {
                //write!(f, "{val:5.3} ");
            }
            //writeln!(f, "");
        }
        Ok(())
    }
}
