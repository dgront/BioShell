use crate::msa::MSA;
use crate::sequence::{count_identical, len_ungapped, ProfileColumnOrder, remove_gaps, Sequence, SequenceProfile};
use crate::SequenceError;

/// Returns the most probable sequence.
///
/// Internally, it uses [`most_probable_sequence()`](SequenceProfile::most_probable_sequence) method
/// from the [`SequenceProfile`] struct. Results in [`SequenceError`] when the [`MSA`] contains no sequences
pub fn most_probable_sequence(msa: &MSA) -> Result<Sequence, SequenceError> {

    if msa.n_seq()==0 { return Err(SequenceError::NoInputSequences) }

    let profile = SequenceProfile::new(ProfileColumnOrder::aa_standard_gapped(), &msa);
    Ok(profile.most_probable_sequence())
}

/// Returns the longest sequence in the given [`MSA`].
///
/// The newly returned object is a copy of the longest sequence with gaps removed. Results in
/// [`SequenceError`] when the [`MSA`] contains no sequences.
///
/// # Example
/// ```
/// # use bioshell_seq::sequence::FastaIterator;
/// # use bioshell_seq::SequenceError;
/// # use std::io::BufReader;
/// # fn main() -> Result<(), SequenceError> {
/// use bioshell_seq::msa::{longest_sequence, MSA};
/// let fasta_input="
/// > syb:TZ53_18285
/// --MSQKLKVVIDKAACCGYGVCAEICPQVYKLDANGIVYVDDEIV-----
/// > 646611275
/// ------MPAKVDENLCTGCGLCEEICPEVFKLDENGISRVVGDCE-----
/// > EGCR1_03845
/// ------MKCEIIPERCIACGLCQTIAPEIFDYTDDGLVLFVGEPEATHEF
/// > SPSE_0307
/// ---MMSYYAYVDRDMCIACSACGAAAPRLFRYDAQGIAYMCLDCNSGTAQ
/// > BBR47_24240
/// MGGETTMTTWVDKDTCIACGACGATAPDVFDYDDEGLAFNKLDDNANSVE
/// ";
/// let mut reader = BufReader::new(fasta_input.as_bytes());
/// let msa = MSA::from_fasta_reader(&mut reader)?;
/// let longest = longest_sequence(&msa)?;
/// assert_eq!(longest.seq(), b"MGGETTMTTWVDKDTCIACGACGATAPDVFDYDDEGLAFNKLDDNANSVE");
/// # Ok(())
/// # }
/// ```
pub fn longest_sequence(msa: &MSA) -> Result<Sequence, SequenceError> {

    if msa.n_seq()==0 { return Err(SequenceError::NoInputSequences) }

    let mut longest_seq = &msa.sequences()[0];
    let mut longest_len = len_ungapped(longest_seq);
    for s in msa.sequences() {
        let l = len_ungapped(s);
        if l> longest_len {
            longest_len = l;
            longest_seq = s;
        }
    }
    let mut seq = longest_seq.clone();
    remove_gaps(&mut seq);

    return Ok(seq);
}

/// Returns the longest sequence in the given [`MSA`].
///
/// The newly returned object is a copy of the longest sequence with gaps removed. Results in
/// [`SequenceError`] when the [`MSA`] contains no sequences.
///
/// # Example
/// ```
/// # use bioshell_seq::sequence::FastaIterator;
/// # use bioshell_seq::SequenceError;
/// # use std::io::BufReader;
/// # fn main() -> Result<(), SequenceError> {
/// use bioshell_seq::msa::{medoid_sequence, MSA};
/// let fasta_input="
/// > syb:TZ53_18285
/// --MSQKLKVVIDKAACCGYGVCAEICPQVYKLDANGIVYVDDEIV-----
/// > 646611275
/// ------MPAKVDENLCTGCGLCEEICPEVFKLDENGISRVVGDCE-----
/// > EGCR1_03845
/// ------MKCEIIPERCIACGLCQTIAPEIFDYTDDGLVLFVGEPEATHEF
/// > SPSE_0307
/// ---MMSYYAYVDRDMCIACSACGAAAPRLFRYDAQGIAYMCLDCNSGTAQ
/// > BBR47_24240
/// MGGETTMTTWVDKDTCIACGACGATAPDVFDYDDEGLAFNKLDDNANSVE
/// ";
/// let mut reader = BufReader::new(fasta_input.as_bytes());
/// let msa = MSA::from_fasta_reader(&mut reader)?;
/// let longest = medoid_sequence(&msa)?;
/// assert_eq!(longest.seq(), b"MPAKVDENLCTGCGLCEEICPEVFKLDENGISRVVGDCE");
/// # Ok(())
/// # }
/// ```
pub fn medoid_sequence(msa: &MSA) -> Result<Sequence, SequenceError> {

    if msa.n_seq()==0 { return Err(SequenceError::NoInputSequences) }

    let mut row_mins = vec![usize::MAX; msa.n_seq()];
    for (i,si) in msa.sequences().iter().enumerate() {
        for (j,sj) in msa.sequences().iter().enumerate() {
            if i == j { break; }
            let ids= count_identical(&sj, si).unwrap();
            eprintln!("{} {} {}", i, j, ids);
            row_mins[i] = row_mins[i].min(ids);
            row_mins[j] = row_mins[j].min(ids);
        }
    }

    let mut best_seq = 0;
    let mut best_val = row_mins[0];

    for i in 1..msa.n_seq() {
        if row_mins[i] > best_val {
            best_val = row_mins[i];
            best_seq = i;
        }
    }

    let mut seq = msa.sequences()[best_seq].clone();
    remove_gaps(&mut seq);

    return Ok(seq);
}