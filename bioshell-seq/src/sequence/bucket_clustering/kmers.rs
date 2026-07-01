use crate::SequenceError;
use crate::SequenceError::InvalidOneLetterCode;
use crate::chemical::standard_letter_to_index;

/// A k-mer encodes a contiguous substring of length k.
///
/// Here we encode k-mers as u32 integers, where each symbol is represented by 5 bits
/// (allowing for up to 32 distinct symbols). This encoding allows to hash k-residue sequence fragments
pub type Kmer = u32;

/// Generate sorted, deduplicated u32-backed k-mers.
///
/// Requirements:
/// - k in 1..=6
/// - every symbol in `seq` must be in 0..=31
/// - different k values should not be mixed unless you encode k separately
pub fn generate_kmers(seq: &[u8], k: usize) -> Result<Vec<u32>, SequenceError> {
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
        let x = standard_letter_to_index(seq[i]).map_err(|_| InvalidOneLetterCode {
            aa_code: seq[i] as char,
            sequence: String::from_utf8_lossy(seq).into_owned(),
        })?;
        assert!(x <= MAX_SYMBOL, "symbol value {x} exceeds maximum allowed value 31");

        code = ((code << BITS_PER_SYMBOL) | x as u32) & mask;
        if i + 1 >= k { kmers.push(code); }
    }

    kmers.sort_unstable();
    kmers.dedup();

    return Ok(kmers);
}


/// Count the number of shared k-mers between two sorted k-mer lists.
///
/// The resulting count can be used to estimate sequence similarity bounds with
/// [`kmer_identity_bounds()`].
#[inline(never)]
pub fn count_intersection_sorted(a: &[u32], b: &[u32]) -> usize {
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

/// Compute bounds on sequence identity based on shared and different k-mers.
///
/// Given two amino acid sequences, we can compute the number of ``different_kmers``
/// that are present in one sequence but not the other. Given the k-mer length and `min_seq_len`
/// - the length of the shorter sequence (used in denominator), we can estimate bounds on the sequence identity.
///
/// The function returns a tuple of (`lower_bound`, `upper_bound`) on the sequence identity fraction.
pub fn kmer_identity_bounds(different_kmers: usize, kmer_len: usize, min_seq_len: usize) -> (f32, f32) {
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
pub fn suggest_word_length(identity_level: f32) -> usize {
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