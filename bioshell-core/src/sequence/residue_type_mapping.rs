/// Provides `u8` (one byte) index for a residue type, e.g. for a standard amino acid.
///
/// [`ResidueTypeOrder`](ResidueTypeOrder) defines the order of columns in a sequence profile, e.g.
/// according to the standard NCBI's order, amino acids are listed alphabetically by their full name,
/// which results in the following order or one-letter codes: `ARNDCQEGHILKMFPSTWYVX`.
///
#[derive(Clone, Debug)]
pub struct ResidueTypeOrder {
    index_to_aa: Vec<u8>,
    aa_to_index: Vec<usize>,
}

impl ResidueTypeOrder {

    /// Creates a new mapping for a given order of letters
    pub fn new(chars_ordered: &str) -> ResidueTypeOrder {
        let index_to_aa = chars_ordered.as_bytes().to_vec();
        let mut aa_to_index: Vec<usize> = vec![0; 128];
        for (i, aai) in chars_ordered.chars().enumerate() {
            aa_to_index[aai as usize] = i;
        }
        ResidueTypeOrder {index_to_aa, aa_to_index}
    }

    /// Creates a new mapping for amino acids in the NCBI's order: `ARNDCQEGHILKMFPSTWYVX`
    pub fn aa_standard() -> ResidueTypeOrder { ResidueTypeOrder::new("ARNDCQEGHILKMFPSTWYVX") }

    /// Creates a new mapping for amino acids and the gap symbol in the NCBI's order: `ARNDCQEGHILKMFPSTWYVX-`
    pub fn aa_standard_gapped() -> ResidueTypeOrder { ResidueTypeOrder::new("ARNDCQEGHILKMFPSTWYVX-") }

    /// Returns the size of this mapping i.e. the number of residue types mapped
    pub fn size(&self) -> usize { self.index_to_aa.len() }

    /// Convert a given amino acid (or nucleotide) character to its order index
    pub fn letter_to_index(&self, aa: &char) -> usize {
        self.aa_to_index[(*aa as u8) as usize]
    }

    /// Convert a given amino acid (or nucleotide) byte to its order index
    pub fn type_to_index(&self, aa: &u8) -> usize {
        self.aa_to_index[*aa as usize]
    }

    /// Converts a residue type order index into its letter
    pub fn index_to_letter(&self, res_id:u8) -> char { self.index_to_aa[res_id as usize] as char }

    /// Converts a residue type order index into its u8 type
    pub fn index_to_type(&self, res_id:u8) -> u8 { self.index_to_aa[res_id as usize] }
}