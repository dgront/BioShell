//! Provides [`Sequence`](crate::sequence::Sequence) and [`SequenceProfile`](crate::sequence::SequenceProfile) stucts
//! as well as numerous utility functions to operate on them


mod sequence;
mod sequence_profile;
mod residue_type_mapping;
mod sequence_filters;
mod sequence_reporter;
pub use sequence_reporter::*;

use log::info;
use bioshell_io::open_file;
pub use sequence::*;
pub use sequence_filters::*;
pub use sequence_profile::*;
pub use residue_type_mapping::*;
use crate::SequenceError;

/// Returns a list of Sequences for a given input string.
///
/// If the input string contains a dot character, it's assumed to be a file name. This function
/// attempts to open that file as FASTA and load all sequences stored.
/// Otherwise, it's assumed the string is an amino acid sequence by itself; it's converted into sequence and returned
/// as the only element of a vector.
pub fn load_sequences(seq_or_fname: &String, seq_name: &str) -> Result<Vec<Sequence>, SequenceError> {
    if seq_or_fname.contains(".") {
        let reader = open_file(&seq_or_fname)?;
        let seq_iter = FastaIterator::new(reader);
        let out: Vec<Sequence> = seq_iter.collect::<Vec<Sequence>>();
        info!("{}",format!("{} sequences loaded from {}", out.len(), &seq_or_fname));
        return Ok(out);
    }

    return Ok(vec![Sequence::from_str(seq_name, seq_or_fname)]);
}

/// Returns the first word of a sequence description string.
///
/// This function is used to extract the ID of a sequence from it's description, assuming it's the first word
///
/// # Examples
/// ```
/// use bioshell_seq::sequence::first_word_of_description;
/// let description = "1ABC A hypothetical protein from E. coli";
/// let id = first_word_of_description(description);
/// assert_eq!(id, "1ABC");
/// let description = "NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome";
/// let id = first_word_of_description(description);
/// assert_eq!(id, "NC_000913.3");
/// ```
pub fn first_word_of_description(description: &str) -> String {
    let v: Vec<&str> = description.split_whitespace().collect();

    return String::from(v[0]);
}