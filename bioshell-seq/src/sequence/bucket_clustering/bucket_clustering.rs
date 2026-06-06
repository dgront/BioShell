
use rand::seq::SliceRandom;
use rayon::prelude::*;
use rand::thread_rng;
use std::time::Instant;
use log::{debug, error, info};
use crate::alignment::{aligned_sequences, GlobalAligner};

use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::{count_identical, Sequence};

use crate::sequence::bucket_clustering::kmers::*;
use crate::SequenceError;

// type Cluster = Vec<usize>;

pub struct BucketClustering {
    pub id_level: f32,
    pub word_size: usize,
    sequences: Vec<Sequence>,
    sequence_order: Vec<usize>,
    kmer_sets: Vec<Vec<Kmer>>,
    aligner: GlobalAligner<SequenceSimilarityScore>,
    scoring: SequenceSimilarityScore,
    n_aligned: usize,
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

pub fn bucket_clustering<'a>(sequences: &'a [Sequence], id_level: f32) -> Result<Vec<Vec<&'a Sequence>>, SequenceError> {

    let mut bc = BucketClustering::new(sequences.to_vec(), id_level)?;
    let clustering = bc.run();
    println!("{:?}", clustering);
    return Ok(clustering.into_iter().map(|cluster| {
                    cluster.members.into_iter().map(|idx| &sequences[idx]).collect()
                }).collect());
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

        // --- prepare for sequence alignment
        let longest_length = sequences.iter().map(|s| s.len()).max().unwrap();
        let aligner = GlobalAligner::new(longest_length);
        let scoring = SequenceSimilarityScore::new(SubstitutionMatrixList::BLOSUM62);

        Ok(Self { id_level, word_size, sequences, sequence_order, kmer_sets, aligner, scoring, n_aligned: 0 })
    }

    pub fn run(&mut self) -> Clustering {

        let singleton_clusters: Vec<Cluster> = self.sequence_order.iter().map(|&i| Cluster::new(i)).collect();
        let mut empty_clustering: Clustering = Vec::new();
        self.merge(&mut empty_clustering, &singleton_clusters);

        return empty_clustering;
        // let one_percent = (self.sequences.len() as f64 * 0.01) as usize;
        // let mut cnt = 0;
    }

    fn merge(&mut self, clusters1: &mut Clustering, clusters2: &Clustering) {

        for b_cluster in clusters2 {
            let mut assigned = false;
            for a_cluster in clusters1.iter_mut() {
                let identity = self.sequence_identity(a_cluster.representative, b_cluster.representative);
                if identity >= self.id_level {
                    a_cluster.extend(b_cluster);
                    assigned = true;
                    break;
                }
            }
            if !assigned {
                clusters1.push(b_cluster.clone());
            }
        }
    }

    fn sequence_identity(&mut self, seq_id1: usize, seq_id2: usize) -> f32 {

        let rep = &self.sequences[seq_id1];
        let candidate = &self.sequences[seq_id2];

        let rep_kmer_set = &self.kmer_sets[seq_id1];
        let candidate_kmers = &self.kmer_sets[seq_id2];

        let n_shared = count_intersection_sorted(candidate_kmers, rep_kmer_set);
        let different = candidate_kmers.len().saturating_sub(n_shared);

        let shorter_len = candidate.len().min(rep.len());
        let (lower, upper) = kmer_identity_bounds(different, self.word_size, shorter_len);

        if lower >= self.id_level { return lower; }
        else if upper < self.id_level {
            return upper;
        } else {
            // Ambiguous: fall back to exact sequence identity check
            debug!("Sequence identity in range {:.2} {:.2} for: {} {}, computing alignment",lower,upper,rep.description_n(10),candidate.description_n(10));

            self.scoring.template_from_sequence(candidate);
            self.scoring.query_from_sequence(rep);

            self.aligner.align(&self.scoring, -11, -1);
            self.n_aligned += 1;

            let path = self.aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, rep, candidate, '-');

            let n_identical = count_identical(&ali_q, &ali_t).unwrap();
            return n_identical as f32 / shorter_len as f32;
        }
    }
}
//
// /// Groups sequences into clusters based on pairwise sequence identity using a
// /// CD-HIT-style greedy clustering algorithm with k-mer-based prefiltering.
// ///
// /// This function partitions sequences into non-overlapping clusters ("buckets"),
// /// where each sequence in a cluster shares a minimum sequence identity with the
// /// cluster's representative sequence. It uses a fast k-mer intersection filter to
// /// avoid unnecessary pairwise alignments, followed by identity estimation using
// /// ungapped residue matching.
// ///
// /// ### Arguments
// ///
// /// * `sequences` - A vector of `Sequence` objects to be clustered. Each sequence
// ///   must implement the `SequenceLike` trait, which provides access to the sequence data as `&[u8]`.
// /// * `id_level` - The minimum pairwise sequence identity (between 0.0 and 1.0) required
// ///   for two sequences to be placed in the same cluster.
// ///
// /// ### Returns
// ///
// /// A vector of clusters, where each cluster is a `Vec` of references to sequences
// /// from the input. Each cluster is represented by its first sequence (the longest
// /// unassigned sequence at insertion time).
// ///
// pub fn bucket_clustering<'a>(sequences: &'a [Sequence], id_level: f32) -> Vec<Vec<&'a Sequence>> {
//     let indexes: Vec<usize> = (0..sequences.len()).collect();
//     let clusters = bucket_clustering_impl(sequences, &indexes, id_level);
//     clusters_to_refs(sequences, clusters)
// }
//
// pub fn bucket_clustering_n<'a>(sequences: &'a [Sequence], id_level: f32, n: usize) -> Vec<Vec<&'a Sequence>>
//                 where Sequence: Sync {
//     let indexes: Vec<usize> = (0..sequences.len()).collect();
//     let clusters = bucket_clustering_n_impl(sequences, indexes, id_level, n);
//     clusters_to_refs(sequences, clusters)
// }
//
// fn bucket_clustering_n_impl(sequences: &[Sequence], mut indexes: Vec<usize>, id_level: f32, n: usize) -> Clustering
//                     where Sequence: Sync {
//
//     if indexes.is_empty() { return Vec::new(); }
//
//     if n <= 1 || indexes.len() <= 1 {
//         return bucket_clustering_impl(sequences, &indexes, id_level);
//     }
//
//     let n = n.min(indexes.len());
//
//     indexes.shuffle(&mut thread_rng());
//
//     let buckets = split_round_robin(indexes, n);
//
//     let local_clusterings: Vec<Clustering> = buckets
//         .into_par_iter()
//         .map(|bucket| bucket_clustering_impl(sequences, &bucket, id_level))
//         .collect();
//
//     let representative_sets: Vec<Vec<usize>> = local_clusterings
//         .iter()
//         .map(representatives_from_clustering)
//         .filter(|reps| !reps.is_empty())
//         .collect();
//
//     let representative_clusters =
//         merge_representative_sets(sequences, representative_sets, id_level);
//
//     repack_local_clusters(
//         sequences.len(),
//         local_clusterings,
//         representative_clusters,
//     )
// }
//
//
// fn repack_local_clusters(n_sequences: usize, local_clusterings: Vec<Clustering>, representative_clusters: Clustering) -> Clustering {
//
//     let mut rep_to_final_cluster: Vec<Option<usize>> = vec![None; n_sequences];
//
//     for (cluster_id, rep_cluster) in representative_clusters.iter().enumerate() {
//         for &rep_idx in rep_cluster {
//             rep_to_final_cluster[rep_idx] = Some(cluster_id);
//         }
//     }
//
//     let mut final_clusters: Clustering = vec![Vec::new(); representative_clusters.len()];
//
//     for local_clustering in local_clusterings {
//         for local_cluster in local_clustering {
//             if let Some(&local_rep) = local_cluster.first() {
//                 let final_cluster_id = rep_to_final_cluster[local_rep]
//                     .expect("missing representative in final clustering");
//
//                 final_clusters[final_cluster_id].extend(local_cluster);
//             }
//         }
//     }
//
//     final_clusters
// }
//
// fn bucket_clustering_impl(sequences: &[Sequence], indexes: &[usize], id_level: f32) -> Clustering {
//
//     // --- nothing to cluster
//     if indexes.is_empty() {
//         return Vec::new();
//     }
//
//     // --- start a timer
//     let start = Instant::now();
//
//     // --- prepare for sequence alignment
//     let longest_length = indexes
//         .iter()
//         .map(|&i| sequences[i].len())
//         .max()
//         .unwrap();
//
//     let mut aligner = GlobalAligner::new(longest_length);
//     let mut scoring = SequenceSimilarityScore::new(SubstitutionMatrixList::BLOSUM62);
//     let mut n_aligned = 0;
//
//     let word_size = suggest_word_length(id_level);
//
//     // --- sort indexes by sequence length in descending order,
//     // --- as in the original CD-HIT implementation
//     let mut sorted_indices: Vec<usize> = indexes.to_vec();
//
//     sorted_indices.sort_by_key(|&i| -(sequences[i].seq().len() as isize));
//
//     // --- precompute k-mer sets for all selected sequences
//     //
//     // kmer_sets is indexed by global sequence index, not by local position.
//     // This keeps the rest of the code simple.
//     let mut kmer_sets: Vec<Option<Vec<Kmer>>> = vec![None; sequences.len()];
//
//     for &i in indexes {
//         match generate_kmers(sequences[i].seq(), word_size) {
//             Ok(kmers) => {
//                 kmer_sets[i] = Some(kmers);
//             }
//             Err(err) => {
//                 error!("{}", err);
//             }
//         }
//     }
//
//     let mut clusters: Clustering = Vec::new();
//     let mut representatives: Vec<usize> = Vec::new();
//     let mut rep_kmers: Vec<&Vec<Kmer>> = Vec::new();
//
//     let one_percent = (sorted_indices.len() as f64 * 0.01) as usize;
//     let mut cnt = 0;
//
//     for &i in &sorted_indices {
//         let candidate = &sequences[i];
//
//         let Some(candidate_kmers) = kmer_sets[i].as_ref() else {
//             // Sequence was skipped because k-mer generation failed.
//             continue;
//         };
//
//         let mut assigned = false;
//
//         for (cluster_idx, (&rep_idx, rep_kmer_set)) in representatives
//             .iter()
//             .zip(&rep_kmers)
//             .enumerate()
//         {
//             let rep = &sequences[rep_idx];
//
//             let n_shared = count_intersection_sorted(candidate_kmers, rep_kmer_set);
//             let different = candidate_kmers.len().saturating_sub(n_shared);
//
//             let shorter_len = candidate.len().min(rep.len());
//             let (lower, upper) = kmer_identity_bounds(different, word_size, shorter_len);
//
//             if lower >= id_level {
//                 // Fast inclusion based on lower bound - they must be in the same cluster
//                 clusters[cluster_idx].push(i);
//                 assigned = true;
//                 break;
//             } else if upper < id_level {
//                 // Fast exclusion based on upper bound - they must be in a different cluster
//                 continue;
//             } else {
//                 debug!(
//                     "Sequence identity in range {:.2} {:.2} for: {} {}, computing alignment",lower,upper,rep.description_n(10),candidate.description_n(10));
//
//                 // Ambiguous: fall back to exact sequence identity check
//                 scoring.template_from_sequence(candidate);
//                 scoring.query_from_sequence(rep);
//
//                 aligner.align(&scoring, -11, -1);
//                 n_aligned += 1;
//
//                 let path = aligner.backtrace();
//                 let (ali_q, ali_t) = aligned_sequences(&path, rep, candidate, '-');
//
//                 let n_identical = count_identical(&ali_q, &ali_t).unwrap();
//                 let identity = n_identical as f32 / shorter_len as f32;
//
//                 if identity >= id_level {
//                     clusters[cluster_idx].push(i);
//                     assigned = true;
//                     break;
//                 }
//             }
//         }
//
//         if !assigned {
//             clusters.push(vec![i]);
//             representatives.push(i);
//             let kmers = kmer_sets[i].as_ref().expect("representative must have k-mers");
//             rep_kmers.push(kmers);
//         }
//
//         if one_percent > 0 {
//             cnt += 1;
//
//             if cnt % one_percent == 1 {
//                 info!("{}% sequences split into buckets",(cnt - 1) / one_percent);
//             }
//         }
//     }
//
//     info!(
//         "{} sequences clustered into {} buckets in {:.3?}, alignment called {} times",
//         indexes.len(),
//         clusters.len(),
//         start.elapsed(),
//         n_aligned
//     );
//
//     clusters
// }
//
// /// Split a list of sequence indices into `n` buckets in a round-robin fashion.
// fn split_round_robin(indexes: Vec<usize>, n: usize) -> Vec<Vec<usize>> {
//
//     let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n];
//     for (i, idx) in indexes.into_iter().enumerate() { buckets[i % n].push(idx); }
//
//     return buckets;
// }
//
// /// Takes a representative sequence from each cluster
// fn representatives_from_clustering(clustering: &Clustering) -> Vec<usize> {
//
//     clustering
//         .iter()
//         .filter_map(|cluster| cluster.first().copied())
//         .collect()
// }
//
//
// fn merge_representative_sets(sequences: &[Sequence], representative_sets: Vec<Vec<usize>>, id_level: f32) -> Clustering
//     where Sequence: Sync {
//     let mut cluster_sets: Vec<Clustering> = representative_sets
//         .into_iter()
//         .map(|set| set.into_iter().map(|idx| vec![idx]).collect())
//         .collect();
//
//     if cluster_sets.is_empty() { return Vec::new(); }
//
//     while cluster_sets.len() > 1 {
//         let current_sets = std::mem::take(&mut cluster_sets);
//
//         let mut pair_jobs: Vec<(Clustering, Clustering)> = Vec::new();
//         let mut carried: Option<Clustering> = None;
//
//         let mut iter = current_sets.into_iter();
//
//         while let Some(left) = iter.next() {
//             if let Some(right) = iter.next() {
//                 pair_jobs.push((left, right));
//             } else {
//                 carried = Some(left);
//             }
//         }
//
//         cluster_sets = pair_jobs
//             .into_par_iter()
//             .map(|(left_clusters, right_clusters)| {
//                 merge_two_representative_clusterings(
//                     sequences,
//                     left_clusters,
//                     right_clusters,
//                     id_level,
//                 )
//             })
//             .collect();
//
//         if let Some(unpaired) = carried {
//             cluster_sets.push(unpaired);
//         }
//     }
//
//     cluster_sets.pop().unwrap()
// }
//
// fn merge_two_representative_clusterings(sequences: &[Sequence], left_clusters: Clustering, right_clusters: Clustering, id_level: f32) -> Clustering {
//     let mut all_clusters = left_clusters;
//     all_clusters.extend(right_clusters);
//
//     let reps: Vec<usize> = all_clusters
//         .iter()
//         .filter_map(|cluster| cluster.first().copied())
//         .collect();
//
//     let rep_clusters = bucket_clustering_impl(sequences, &reps, id_level);
//
//     let mut rep_to_old_cluster: Vec<Option<usize>> = vec![None; sequences.len()];
//
//     for (old_cluster_id, old_cluster) in all_clusters.iter().enumerate() {
//         if let Some(&rep) = old_cluster.first() {
//             rep_to_old_cluster[rep] = Some(old_cluster_id);
//         }
//     }
//
//     let mut merged_clusters: Clustering = Vec::with_capacity(rep_clusters.len());
//
//     for rep_cluster in rep_clusters {
//         let mut merged_cluster = Vec::new();
//
//         for rep in rep_cluster {
//             let old_cluster_id = rep_to_old_cluster[rep]
//                 .expect("representative not found during representative merge");
//
//             merged_cluster.extend_from_slice(&all_clusters[old_cluster_id]);
//         }
//
//         merged_clusters.push(merged_cluster);
//     }
//
//     merged_clusters
// }
//
// /// Convert clusters of sequence indices back to clusters of sequence references.
// ///
// /// This back-mapping happens after the clustering is finished
// fn clusters_to_refs<'a>(sequences: &'a [Sequence], clusters: Clustering) -> Vec<Vec<&'a Sequence>> {
//     clusters.into_iter().map(|cluster| {
//             cluster.into_iter().map(|idx| &sequences[idx]).collect()
//         }).collect()
// }
//
