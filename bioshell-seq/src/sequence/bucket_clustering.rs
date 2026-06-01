use log::{debug, error, info};
use crate::alignment::{aligned_sequences, GlobalAligner};
use crate::chemical::standard_letter_to_index;
use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::{count_identical, Sequence};
use crate::SequenceError;
use crate::SequenceError::InvalidOneLetterCode;

type Kmer = u32;


/// Generate sorted, deduplicated u32-backed k-mers.
///
/// Requirements:
/// - k in 1..=6
/// - every symbol in `seq` must be in 0..=31
/// - different k values should not be mixed unless you encode k separately
fn generate_kmers(seq: &[u8], k: usize) -> Result<Vec<u32>, SequenceError> {
    const BITS_PER_SYMBOL: usize = 5;
    const MAX_SYMBOL: u8 = 31;
    const MAX_K_U32: usize = 6;

    if k == 0 || k > MAX_K_U32 || seq.len() < k {
        return Ok(Vec::new());
    }

    let mask: u32 = (1u32 << (BITS_PER_SYMBOL * k)) - 1;
    let mut code = 0u32;

    let mut kmers = Vec::with_capacity(seq.len() - k + 1);

    for i in 0..seq.len() {
        let x = standard_letter_to_index(seq[i])?;

        assert!(x <= MAX_SYMBOL, "symbol value {x} exceeds maximum allowed value 31");

        code = ((code << BITS_PER_SYMBOL) | x as u32) & mask;
        if i + 1 >= k { kmers.push(code); }
    }

    kmers.sort_unstable();
    kmers.dedup();

    return Ok(kmers);
}

#[inline(never)]
fn count_intersection_sorted(a: &[u32], b: &[u32]) -> usize {
    let mut i = 0;
    let mut j = 0;
    let mut count = 0;

    let a_len = a.len();
    let b_len = b.len();

    while i < a_len && j < b_len {
        let av = unsafe { *a.get_unchecked(i) };
        let bv = unsafe { *b.get_unchecked(j) };

        if av < bv {
            i += 1;
        } else if av > bv {
            j += 1;
        } else {
            count += 1;
            i += 1;
            j += 1;
        }
    }

    return count;
}

/// Groups sequences into clusters based on pairwise sequence identity using a
/// CD-HIT-style greedy clustering algorithm with k-mer-based prefiltering.
///
/// This function partitions sequences into non-overlapping clusters ("buckets"),
/// where each sequence in a cluster shares a minimum sequence identity with the
/// cluster's representative sequence. It uses a fast k-mer intersection filter to
/// avoid unnecessary pairwise alignments, followed by identity estimation using
/// ungapped residue matching.
///
/// ### Arguments
///
/// * `sequences` - A vector of `Sequence` objects to be clustered. Each sequence
///   must implement the `SequenceLike` trait, which provides access to the sequence data as `&[u8]`.
/// * `id_level` - The minimum pairwise sequence identity (between 0.0 and 1.0) required
///   for two sequences to be placed in the same cluster.
///
/// ### Returns
///
/// A vector of clusters, where each cluster is a `Vec` of references to sequences
/// from the input. Each cluster is represented by its first sequence (the longest
/// unassigned sequence at insertion time).
///
pub fn bucket_clustering<'a>(sequences: &'a Vec<Sequence>, id_level: f32) -> Vec<Vec<&'a Sequence>> {

    // --- prepare for sequence alignment
    let longest_length = sequences.iter().map(|s|s.len()).max().unwrap();
    let mut aligner = GlobalAligner::new(longest_length);
    let mut scoring = SequenceSimilarityScore::new( SubstitutionMatrixList::BLOSUM62);
    let mut n_aligned = 0;

    let word_size = suggest_word_length(id_level);

    // --- sort sequences by length in descending order, as in the original CD-HIT implementation
    let mut sorted_indices: Vec<usize> = (0..sequences.len()).collect();
    sorted_indices.sort_by_key(|&i| -(sequences[i].seq().len() as isize));

    // Precompute k-mer sets for all sequences
    let mut kmer_sets: Vec<Vec<Kmer>> = Vec::with_capacity(sequences.len());
    for s in sequences {
        // --- here we handle the case when generate_kmers() found illegal residue code
        match generate_kmers(s.seq(), word_size) {
            Ok(kmers) => kmer_sets.push(kmers),
            Err(InvalidOneLetterCode { aa_code, .. }) => {
                let new_err = InvalidOneLetterCode { aa_code, sequence: String::from_utf8_lossy(s.seq()).into_owned() };
                error!("{}", new_err);
            }
            _ => {}
        }
    }

    let mut clusters: Vec<Vec<&'a Sequence>> = Vec::new();
    let mut representatives: Vec<&'a Sequence> = Vec::new();
    let mut rep_kmers: Vec<&Vec<Kmer>> = Vec::new();
    let one_percent = (sorted_indices.len() as f64 * 0.01) as usize;
    let mut cnt = 0;
    for &i in &sorted_indices {
        let candidate = &sequences[i];
        let candidate_kmers = &kmer_sets[i];
        let mut assigned = false;

        for (cluster_idx, (rep, rep_kmer_set)) in representatives.iter().zip(&rep_kmers).enumerate() {
            let n_shared = count_intersection_sorted(candidate_kmers, rep_kmer_set);
            let different = candidate_kmers.len().saturating_sub(n_shared);

            let shorter_len = candidate.len().min(rep.len());
            let (lower, upper) = kmer_identity_bounds(different, word_size, shorter_len);
            if lower >= id_level {
                // Fast inclusion based on lower bound - they must be in the same cluster
                clusters[cluster_idx].push(candidate);
                assigned = true;
                break;
            } else if upper < id_level {
                // Fast exclusion based on upper bound - they must be in a different cluster
                continue;
            } else {
                debug!("Sequence identity in range {:.2} {:.2} for: {} {}, computing alignment", lower, upper, rep.description_n(10), candidate.description_n(10));
                // Ambiguous: fall back to exact sequence identity check
                scoring.template_from_sequence(candidate);
                scoring.query_from_sequence(rep);
                aligner.align(&scoring, -11, -1);
                n_aligned += 1;
                let path = aligner.backtrace();
                let (ali_q, ali_t) = aligned_sequences(&path, rep, candidate, '-');
                let n_identical = count_identical(&ali_q, &ali_t).unwrap();
                let identity = n_identical as f32 / shorter_len as f32;
                if identity >= id_level {
                    clusters[cluster_idx].push(candidate);
                    assigned = true;
                    break;
                }
            }
        }

        if !assigned {
            clusters.push(vec![candidate]);
            representatives.push(candidate);
            rep_kmers.push(&kmer_sets[i]);
        }
        if one_percent > 0 {
            cnt += 1;
            if cnt % one_percent == 1 {
                info!("{}% sequences split into buckets", (cnt-1) / one_percent);
            }
        }
    }
    debug!("Sequence  alignment called {} times", n_aligned);

    clusters
}

/// Compute bounds on sequence identity based on shared and different k-mers.
/// Assumes that `seq_len` is the length of the shorter sequence (used in denominator).
fn kmer_identity_bounds(different_kmers: usize, kmer_len: usize, min_seq_len: usize) -> (f32, f32) {
    if min_seq_len == 0 { return (0.0, 0.0); }

    // upper bound: best case = each mutation causes k k-mers
    let min_mutations = (different_kmers) / kmer_len + 1; // ceil(d / k)
    let upper_identity = (min_seq_len - min_mutations) as f32 / min_seq_len as f32;

    // lower bound
    let max_mutations = different_kmers + kmer_len - 1;
    let lower_identity = (min_seq_len - max_mutations) as f32 / min_seq_len as f32;

    (lower_identity.max(0.0), upper_identity.min(1.0))
}


/// Suggests a "good" k-mer (word) length for a given identity level.
/// Based on heuristics inspired by CD-HIT.
///
/// # Arguments
/// * `identity_level` - Expected sequence identity (e.g., 0.9 for 90%)
///
/// # Returns
/// Recommended word size (`k`)
fn suggest_word_length(identity_level: f32) -> usize {
    match identity_level {
        x if x >= 0.95 => 6,
        x if x >= 0.90 => 5,
        x if x >= 0.85 => 5,
        x if x >= 0.80 => 4,
        x if x >= 0.75 => 4,
        x if x >= 0.70 => 3,
        x if x >= 0.60 => 3,
        x if x >= 0.50 => 2,
        _ => 1, // fallback for very low identity, but likely not meaningful
    }
}