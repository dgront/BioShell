use bioshell_seq::alignment::{AlignmentReporter, PrintAsPairwise};
use bioshell_seq::sequence::Sequence;

#[test]
fn report_pairwise() {
    let query = Sequence::from_str("query", "AL-IV");
    let template = Sequence::from_str("template", "ALRIV");
    let mut reporter = PrintAsPairwise::new(10);
    reporter.report(&query, &template);
}