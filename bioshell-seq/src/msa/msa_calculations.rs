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
/// [`SequenceError`] when the [`MSA`] contains no sequences
pub fn longest_sequence(msa: &MSA) -> Result<Sequence, SequenceError> {

    if msa.n_seq()==0 { return Err(SequenceError::NoInputSequences) }

    let mut longest_seq = &msa.sequences()[0];
    let mut longest_len = longest_seq.len();
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

pub fn medoid_sequence(msa: &MSA) -> Result<Sequence, SequenceError> {

    if msa.n_seq()==0 { return Err(SequenceError::NoInputSequences) }

    let mut row_mins = vec![usize::MAX; msa.n_seq()];
    for (i,si) in msa.sequences().iter().enumerate() {
        for (j,sj) in msa.sequences().iter().enumerate() {
            if i == j { break; }
            let ids= count_identical(&sj, si).unwrap();
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