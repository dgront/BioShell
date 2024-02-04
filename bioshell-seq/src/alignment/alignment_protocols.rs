use crate::alignment::{aligned_sequences, AlignmentReporter, GlobalAligner};
use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::Sequence;

pub fn align_all_pairs<R: AlignmentReporter>(queries: &Vec<Sequence>, templates: &Vec<Sequence>,
                 matrix: SubstitutionMatrixList, gap_open: i32, gap_extend: i32, lower_triangle: bool,
                 reporter: &mut R) {

    let mut max_length = queries.iter().map(|s| s.len()).max().unwrap();
    max_length = max_length.max(templates.iter().map(|s| s.len()).max().unwrap());
    let mut aligner = GlobalAligner::new(max_length);
    let mut scoring = SequenceSimilarityScore::new( matrix);

    for template in templates {
        scoring.template_from_sequence(template);
        for query in queries {
            if lower_triangle && template==query { break }
            scoring.query_from_sequence(query);
            aligner.align(&scoring, gap_open, gap_extend);
            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, &query, &template, '-');
            reporter.report(&ali_q, &ali_t);
        }
    }
}