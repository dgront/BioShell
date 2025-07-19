use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use log::{debug};
use crate::alignment::{aligned_sequences, GlobalAligner};
use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::{count_identical, Sequence};

/// Wrapper to enable Hash and Eq on &[u8]
#[derive(Debug, Copy, Clone)]
struct Kmer<'a>(&'a [u8]);

impl<'a> PartialEq for Kmer<'a> {
    fn eq(&self, other: &Self) -> bool { self.0 == other.0 }
}

impl<'a> Eq for Kmer<'a> {}

impl<'a> Hash for Kmer<'a> {
    fn hash<H: Hasher>(&self, state: &mut H) { state.write(self.0); }
}

/// Generate k-mers as &[u8] slices wrapped in Kmer<'a>
fn generate_kmers<'a>(seq: &'a [u8], k: usize) -> HashSet<Kmer<'a>> {
    let mut kmers = HashSet::new();
    if seq.len() < k {
        return kmers;
    }
    for i in 0..=seq.len() - k {
        kmers.insert(Kmer(&seq[i..i + k]));
    }
    kmers
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
    let kmer_sets: Vec<HashSet<Kmer<'_>>> = sequences
        .iter()
        .map(|s| generate_kmers(s.seq(), word_size))
        .collect();

    let mut clusters: Vec<Vec<&'a Sequence>> = Vec::new();
    let mut representatives: Vec<&'a Sequence> = Vec::new();
    let mut rep_kmers: Vec<&HashSet<Kmer<'_>>> = Vec::new();

    for &i in &sorted_indices {
        let candidate = &sequences[i];
        let candidate_kmers = &kmer_sets[i];
        let mut assigned = false;

        for (cluster_idx, (rep, rep_kmer_set)) in representatives.iter().zip(&rep_kmers).enumerate() {
            let shared = candidate_kmers
                .intersection(rep_kmer_set)
                .count();
            let different = candidate_kmers.len().saturating_sub(shared);

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