use std::time::Instant;
use log::{debug, info};
use crate::alignment::{aligned_sequences, AlignmentReporter, GlobalAligner};
use crate::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use crate::sequence::Sequence;

pub fn align_all_pairs(queries: &Vec<Sequence>, templates: &Vec<Sequence>,
                 matrix: SubstitutionMatrixList, gap_open: i32, gap_extend: i32, lower_triangle: bool,
                 reporters: &mut Vec<Box<dyn AlignmentReporter>>) {

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
            if lower_triangle && template==query { break }
            scoring.query_from_sequence(query);
            aligner.align(&scoring, gap_open, gap_extend);
            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, &query, &template, '-');
            for reporter in &mut *reporters {
                reporter.report(&ali_q, &ali_t);
                if lower_triangle { reporter.report(&ali_t, &ali_q); }
            }
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