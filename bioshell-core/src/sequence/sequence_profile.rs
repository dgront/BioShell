use std::fmt;
use core::fmt::{Formatter, Display};
use std::fs::File;
use std::collections::HashMap;

use crate::sequence::Sequence;

/// Provides `u8` (one byte) index for a residue type letter, e.g. for a standard amino acid
#[derive(Clone, Debug)]
pub struct ResidueTypeOrder {
    index_to_aa: Vec<char>,
    aa_to_index_map: HashMap<char, u8>,
}

impl ResidueTypeOrder {

    /// Creates a new mapping for a given order of letters
    pub fn new(chars_ordered: &str) -> ResidueTypeOrder {
        let index_to_aa = chars_ordered.chars().collect();
        let mut aa_to_index_map: HashMap<char, u8> = HashMap::new();
        for (i, aai) in chars_ordered.chars().enumerate() {
            aa_to_index_map.insert(aai, i as u8);
        }
        ResidueTypeOrder{index_to_aa, aa_to_index_map}
    }

    /// Creates a new mapping for amino acids in the NCBI's order: `ARNDCQEGHILKMFPSTWYVX`
    pub fn aa_standard() -> ResidueTypeOrder { ResidueTypeOrder::new("ARNDCQEGHILKMFPSTWYVX") }

    /// Returns the size of this mapping i.e. the number of residue types mapped
    pub fn size(&self) -> usize { self.index_to_aa.len() }

    /// Convert a given amino acid (or nucleotide) character to its order index
    pub fn encode_letter(&self, aa: &char) -> u8 {
        let aa_id: &u8 = match self.aa_to_index_map.get(&aa) {
            Some(i) => { i },
            None => {
                eprintln!("Unknown amino acid symbol {}, converted to gap", &aa);
                &self.aa_to_index_map[&'-']
            }
        };

        *aa_id
    }

    /// Converts a residue type order index into its letter
    pub fn decode_letter(&self, res_id:u8) -> char { self.index_to_aa[res_id as usize] }
}

#[derive(Clone, Debug)]
/// Sequence profile describes preference for each of the 20 standard amino acid residue types
/// (possibly also with the gap symbol) at each of the residue positions in the query sequence.
///
/// The preferences are represented as probabilities that add to `1.0` for each sequence position
///
pub struct SequenceProfile {
    mapping: ResidueTypeOrder,
    data: Vec<Vec<f32>>,
}

impl SequenceProfile {
    pub fn new(mapping: ResidueTypeOrder, msa: &Vec<Sequence>) -> SequenceProfile {

        let n = msa[0].len();
        let k = mapping.size();
        let mut data: Vec<Vec<f32>> = vec![vec![0.0; k]; n];

        for sequence in msa {
            if sequence.len() != n {
                eprintln!("\nSequence of incorrect length! Is: {}, should be: {}. The sequence skipped::\n {}\n",
                          sequence.len(), n, sequence);
                continue;
            }            for idx in 0..n {
                let aa_idx = mapping.encode_letter(&sequence.char(idx)) as usize;
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

    pub fn fraction(&self, pos: usize, aa: usize) -> f32 { self.data[pos][aa] }

    pub fn residue_order(&self) -> &ResidueTypeOrder { &self.mapping }

    pub fn len(&self) -> usize { self.data.len() }

}

impl Display for SequenceProfile {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (i, row) in self.data.iter().enumerate() {
            write!(f, "{i:4} ");
            for val in row {
                write!(f, "{val:5.3} ");
            }
            writeln!(f, "");
        }
        Ok(())
    }
}