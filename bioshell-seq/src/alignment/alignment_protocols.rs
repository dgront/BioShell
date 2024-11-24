use std::time::Instant;
use log::{debug, info};
use crate::alignment::{aligned_sequences, AlignmentReporter, GlobalAligner};
use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::Sequence;

/// High-level procedure for aligning all pairs of sequences.
///
/// The procedure aligns each sequence ``i`` from the query set against each sequence ``j`` from the template set,
/// when ``if_triangle_only`` is set to ``true``. When ``if_triangle_only`` is set to ``false``, only ``i < j`` pairs are aligned.
/// Produced alignments are reported to the provided [`reporter`](`AlignmentReporter``).
///
/// # Examples
///
/// To align sequences, first load them into a ``Vec``, e.g. using the [`FastaIterator`](`crate::sequence::FastaIterator``):
/// ```
/// use bioshell_seq::sequence::{Sequence, FastaIterator};
/// let sequences_str: &str =
/// "> 1clf:A
/// AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE
/// > 1dur:A
/// AYVINDSCIACGACKPECPVNCIQEGSIYAIDADSCIDCGSCASVCPVGAPNPED
/// > 1fca:A
/// AYVINEACISCGACEPECPVDAISQGGSRYVIDADTCIDCGACAGVCPVDAPVQA";
///
/// let fasta_iter = FastaIterator::new(sequences_str.as_bytes());
/// let sequences: Vec<Sequence> = fasta_iter.collect();
/// # assert_eq!(sequences.len(), 3);
/// ```
/// ## to align all pairs of sequences from the ``sequences`` vector:
///
/// pass the vector *twice* to the ``align_all_pairs`` function and set the `if_triangle_only` flag to ``false``:
/// ```
/// use bioshell_seq::alignment::{align_all_pairs, SimilarityReport};
/// use bioshell_seq::scoring::SubstitutionMatrixList;
/// # use bioshell_seq::sequence::{Sequence, FastaIterator};
/// # let sequences_str: &str = "> 1clf:A
/// # AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE
/// # > 1dur:A
/// # AYVINDSCIACGACKPECPVNCIQEGSIYAIDADSCIDCGSCASVCPVGAPNPED
/// # > 1fca:A
/// # AYVINEACISCGACEPECPVDAISQGGSRYVIDADTCIDCGACAGVCPVDAPVQA";
/// # let fasta_iter = FastaIterator::new(sequences_str.as_bytes());
/// # let sequences: Vec<Sequence> = fasta_iter.collect();
///
/// let mut reporter = SimilarityReport::new(0);
/// align_all_pairs(&sequences, &sequences, SubstitutionMatrixList::BLOSUM62, -10, -1, true, &mut reporter);
/// ```
///
/// ## to align every sequence from the ``queries`` vector with each sequence from the ``templates`` vector:
///
/// pass the two vectors and set the `if_triangle_only` flag to ``false``:
/// ```
/// use bioshell_seq::alignment::{align_all_pairs, SimilarityReport};
/// use bioshell_seq::scoring::SubstitutionMatrixList;
/// # use bioshell_seq::sequence::{Sequence, FastaIterator};
/// # let sequences_str: &str = "> 1clf:A
/// # AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE
/// # > 1dur:A
/// # AYVINDSCIACGACKPECPVNCIQEGSIYAIDADSCIDCGSCASVCPVGAPNPED
/// # > 1fca:A
/// # AYVINEACISCGACEPECPVDAISQGGSRYVIDADTCIDCGACAGVCPVDAPVQA";
/// # let fasta_iter = FastaIterator::new(sequences_str.as_bytes());
/// # let queries: Vec<Sequence> = fasta_iter.collect();
/// # let fasta_iter = FastaIterator::new(sequences_str.as_bytes());
/// # let templates: Vec<Sequence> = fasta_iter.collect();
///
/// let mut reporter = SimilarityReport::new(0);
/// align_all_pairs(&queries, &templates, SubstitutionMatrixList::BLOSUM62, -10, -1, false, &mut reporter);
/// ```
pub fn align_all_pairs<R: AlignmentReporter>(queries: &Vec<Sequence>, templates: &Vec<Sequence>,
        matrix: SubstitutionMatrixList, gap_open: i32, gap_extend: i32, if_triangle_only: bool, reporter: &mut R) {

    let mut max_length = queries.iter().map(|s| s.len()).max().unwrap();
    max_length = max_length.max(templates.iter().map(|s| s.len()).max().unwrap());
    let mut aligner = GlobalAligner::new(max_length);
    let mut scoring = SequenceSimilarityScore::new( matrix);

    let mut gcups = 0.0;
    let mut n_pairs = 0;
    let start = Instant::now();
    for template in templates {
        scoring.template_from_sequence(template);
        for query in queries {
            if if_triangle_only && template==query { break }
            scoring.query_from_sequence(query);
            aligner.align(&scoring, gap_open, gap_extend);
            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, &query, &template, '-');
            reporter.report(&ali_q, &ali_t);
            if if_triangle_only { reporter.report(&ali_t, &ali_q); }
            gcups += (query.len() * template.len()) as f64;
            n_pairs += 1;
            if n_pairs % 100 == 0 {
                debug!("{}", format!("{} sequence pairs aligned", n_pairs));
            }
        }
    }
    info!("{}", format!("{} sequence pairs aligned in {:?}", n_pairs, start.elapsed()));
    gcups /= 1.0e6;
    gcups /= start.elapsed().as_millis() as f64;
    info!("{}", format!("{} GCUPS", gcups));
}