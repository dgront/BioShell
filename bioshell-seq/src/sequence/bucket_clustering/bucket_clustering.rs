use std::fmt;
use rayon::prelude::*;
use std::time::Instant;
use log::{debug, error, info};
use crate::alignment::{aligned_sequences, GlobalAligner};

use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::{count_identical, Sequence};

use crate::sequence::bucket_clustering::kmers::*;
use crate::SequenceError;

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
/// * `sequences` - A vector of `Sequence` objects to be clustered.
/// * `id_level` - The minimum pairwise sequence identity (between 0.0 and 1.0) required
///   for two sequences to be placed in the same cluster.
///
/// ### Returns
///
/// A vector of clusters, where each cluster is a `Vec` of references to sequences
/// from the input. Each cluster is represented by its first sequence (the longest) sequence.
///
pub fn bucket_clustering<'a>(sequences: &'a [Sequence], id_level: f32) -> Result<Vec<Vec<&'a Sequence>>, SequenceError> {

    let mut bc = BucketClustering::new(sequences.to_vec(), id_level)?;
    let clustering = bc.run();
    return Ok(clustering.into_iter().map(|cluster| {
        cluster.members.into_iter().map(|idx| &sequences[idx]).collect()
    }).collect());
}

/// Parallel version of `bucket_clustering` that uses `n_threads` threads to perform clustering.
///
/// This function uses `n_threads` threads to partition sequences into non-overlapping clusters ("buckets"),
/// where each sequence in a cluster shares a minimum sequence identity with the
/// cluster's representative sequence. It uses a fast k-mer intersection filter to
/// avoid unnecessary pairwise alignments, followed by identity estimation using
/// ungapped residue matching.
///
/// ### Arguments
///
/// * `sequences` - A vector of `Sequence` objects to be clustered.
/// * `id_level` - The minimum pairwise sequence identity (between 0.0 and 1.0) required
///   for two sequences to be placed in the same cluster.
///
/// ### Returns
///
/// A vector of clusters, where each cluster is a `Vec` of references to sequences
/// from the input. Each cluster is represented by its first sequence (the longest) sequence.
pub fn bucket_clustering_n<'a>(sequences: &'a [Sequence], id_level: f32, n_threads: usize) -> Result<Vec<Vec<&'a Sequence>>, SequenceError> {

    let mut bc = BucketClustering::new(sequences.to_vec(), id_level)?;
    let clustering = bc.run_n(n_threads);
    return Ok(clustering.into_iter().map(|cluster| {
        cluster.members.into_iter().map(|idx| &sequences[idx]).collect()
    }).collect());
}

struct BucketClustering {
    pub id_level: f32,
    pub word_size: usize,
    sequences: Vec<Sequence>,
    sequence_order: Vec<usize>,
    kmer_sets: Vec<Vec<Kmer>>,
}

#[derive(Clone, Debug)]
struct Cluster {
    representative: usize,
    members: Vec<usize>,
}

impl Cluster {
    pub fn new(representative: usize) -> Self {
        let members = vec![representative];
        Self { representative, members }
    }
    pub fn extend(&mut self, other: &Cluster) {
        self.members.extend_from_slice(&other.members);
    }
     pub fn push(&mut self, member: usize) {
        self.members.push(member);
    }
}

type Clustering = Vec<Cluster>;

enum SequenceIdentityResult {
    AboveThreshold,
    BelowThreshold,
    Aligned(f32)
}

#[derive(Debug, Clone, Copy, Default)]
struct SequenceIdentityStats {
    pub above_threshold: usize,
    pub below_threshold: usize,
    pub aligned: usize,
}

impl SequenceIdentityStats {
    pub fn add_result(&mut self, result: SequenceIdentityResult) {
        match result {
            SequenceIdentityResult::AboveThreshold => {
                self.above_threshold += 1;
            }
            SequenceIdentityResult::BelowThreshold => {
                self.below_threshold += 1;
            }
            SequenceIdentityResult::Aligned(_) => {
                self.aligned += 1;
            }
        }
    }

    pub fn merge(&mut self, other: SequenceIdentityStats) {
        self.above_threshold += other.above_threshold;
        self.below_threshold += other.below_threshold;
        self.aligned += other.aligned;
    }
}

impl fmt::Display for SequenceIdentityStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let total = self.above_threshold + self.below_threshold + self.aligned;

        write!(f, "above_threshold={}, below_threshold={}, aligned={}, total={}",
            self.above_threshold, self.below_threshold, self.aligned, total)
    }
}

impl BucketClustering {
    pub fn new(sequences: Vec<Sequence>, id_level: f32) -> Result<Self, SequenceError> {

        // --- word size based on the identity level
        let word_size = suggest_word_length(id_level);

        // --- sort indexes by sequence length in descending order, as in the original CD-HIT implementation
        let mut sequence_order: Vec<usize> = (0..sequences.len()).collect();
        sequence_order.sort_by_key(|&i| -(sequences[i].seq().len() as isize));

        // --- precompute k-mer sets for all selected sequences
        //
        // kmer_sets is indexed by global sequence index, not by local position.
        // This keeps the rest of the code simple.
        let mut kmer_sets: Vec<Vec<Kmer>> = Vec::with_capacity(sequences.len());
        for i in 0..sequences.len() {
            kmer_sets.push( generate_kmers(sequences[i].seq(), word_size)? );
        }

        Ok(Self { id_level, word_size, sequences, sequence_order, kmer_sets})
    }

    pub fn run(&mut self) -> Clustering {

        let singleton_clusters: Vec<Cluster> = self.sequence_order.iter().map(|&i| Cluster::new(i)).collect();
        let empty_clustering: Clustering = Vec::new();
        let clusters = self.merge(empty_clustering, singleton_clusters);

        return clusters;
    }

    pub fn run_n(&mut self, n_threads: usize) -> Clustering {

        if n_threads==1 { return self.run(); }

        let singleton_clusters: Vec<Cluster> = self.sequence_order.iter().map(|&i| Cluster::new(i)).collect();

        // sanitize n_threads, calculate chunk size
        let n_threads = n_threads.max(1).min(singleton_clusters.len());
        let chunk_size = singleton_clusters.len().div_ceil(n_threads);

        // Split singleton clusters into chunks and cluster each chunk independently in parallel.
        let mut clusterings: Vec<Clustering> = singleton_clusters
            .par_chunks(chunk_size)
            .map(|chunk| {
                let local_input: Clustering = chunk.to_vec();
                self.merge(Clustering::new(), local_input)
            })
            .collect();

        while clusterings.len() > 1 {
            clusterings = clusterings
                .par_chunks(2)
                .map(|pair| {
                    if pair.len() == 2 {
                        self.merge(pair[0].clone(), pair[1].clone())
                    } else {
                        pair[0].clone()
                    }
                })
                .collect();
        }

        return clusterings.pop().unwrap();
    }

    /// Merge two clusterings by trying to assign each cluster in `clusters2` to a cluster in `clusters1` based on sequence identity criterion
    ///
    /// Consumes both clustering, returns a new Vec<Cluster>.
    fn merge(&self, mut clusters1: Clustering, clusters2: Clustering) -> Clustering {

        info!("merging clusterings of {} and {} clusters", clusters1.len(), clusters2.len());
        let start = Instant::now();

        // --- prepare for sequence alignment
        let mut longest_length: usize = 0;
        for cluster in &clusters1 {
            longest_length = longest_length.max(self.sequences[cluster.representative].len());
        }
        for cluster in &clusters2 {
            longest_length = longest_length.max(self.sequences[cluster.representative].len());
        }
        let mut aligner = GlobalAligner::new(longest_length);
        let mut scoring = SequenceSimilarityScore::new(SubstitutionMatrixList::BLOSUM62);

        // --- reporting progress
        let one_percent = (clusters2.len() as f64 * 0.01) as usize;
        let mut cnt = 0;
        let mut stats = SequenceIdentityStats::default();
        let n_jobs = clusters2.len();
        // --- main loop: for each cluster in the second clustering, try to assign it to a cluster in the first clustering
        for b_cluster in clusters2 {
            let mut assigned = false;
            for a_cluster in clusters1.iter_mut() {
                let identity = self.sequence_identity(&mut scoring, &mut aligner, a_cluster.representative, b_cluster.representative);
                match identity {
                    SequenceIdentityResult::AboveThreshold => {
                        stats.add_result(identity);
                        a_cluster.extend(&b_cluster);
                        assigned = true;
                        break;
                    }
                    SequenceIdentityResult::BelowThreshold => {
                        stats.add_result(identity);
                        continue;
                    }
                    SequenceIdentityResult::Aligned(id) => {
                        stats.add_result(SequenceIdentityResult::Aligned(id));
                        if id >= self.id_level {
                            a_cluster.extend(&b_cluster);
                            assigned = true;
                            break;
                        }
                    }
                }
            }
            if !assigned { clusters1.push(b_cluster); }

            if one_percent > 0 {
                cnt += 1;

                if cnt % 100 == 0 {
                    info!("{} out of {} sequences split into buckets", cnt, n_jobs);
                }
            }
        }
        info!("{} sequences merged into {} buckets in {:.3?}", n_jobs, clusters1.len(), start.elapsed());
        info!("sequence identity stats: {}", stats);

        return clusters1;
    }

    fn sequence_identity(&self, scoring: &mut SequenceSimilarityScore,
            aligner: &mut GlobalAligner<SequenceSimilarityScore>,
                         seq_id1: usize, seq_id2: usize) -> SequenceIdentityResult {

        let rep = &self.sequences[seq_id1];
        let candidate = &self.sequences[seq_id2];

        let rep_kmer_set = &self.kmer_sets[seq_id1];
        let candidate_kmers = &self.kmer_sets[seq_id2];

        let n_shared = count_intersection_sorted(candidate_kmers, rep_kmer_set);
        let different = candidate_kmers.len().saturating_sub(n_shared);

        let shorter_len = candidate.len().min(rep.len());
        let (lower, upper) = kmer_identity_bounds(different, self.word_size, shorter_len);

        if lower >= self.id_level {
            return SequenceIdentityResult::AboveThreshold;
        }
        else if upper < self.id_level {
            return SequenceIdentityResult::BelowThreshold;
        } else {
            // Ambiguous: fall back to exact sequence identity check
            debug!("Sequence identity in range {:.2} {:.2} for: {} {}, computing alignment",lower,upper,rep.description_n(10),candidate.description_n(10));

            scoring.template_from_sequence(candidate);
            scoring.query_from_sequence(rep);

            aligner.align(scoring, -11, -1);
            // self.n_aligned += 1;

            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, rep, candidate, '-');

            let n_identical = count_identical(&ali_q, &ali_t).unwrap();
            return SequenceIdentityResult::Aligned(n_identical as f32 / shorter_len as f32);
        }
    }
}
