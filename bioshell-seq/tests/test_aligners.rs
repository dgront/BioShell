use bioshell_seq::alignment::{GlobalAligner, aligned_strings};
use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};

struct GlobalAlignmentTestCase {
    query: &'static str,
    template: &'static str,
    aligned_query: &'static str,
    aligned_template: &'static str,
    alignment: &'static str,
    score: i32
}

static GLOBAL_CASES: [GlobalAlignmentTestCase; 3] = [
    GlobalAlignmentTestCase {
        query: "A", template: "AW",
        aligned_query: "A-", aligned_template: "AW",
        alignment: "*-", score: -6,
    },
    GlobalAlignmentTestCase {
        query: "AR", template: "ARK",
        aligned_query: "AR-", aligned_template: "ARK",
        alignment: "**-", score: -1,
    },
    GlobalAlignmentTestCase {
        query: "MAVRLLKTHL", template: "MKNITCYL",
        aligned_query: "MAVRLLKTHL", aligned_template: "M--KNITCYL",
        alignment: "*||*******", score: -2,
    },
];

#[test]
fn test_global_aligner() {

    let max_seq_len = 10;
    let mut aligner = GlobalAligner::new(max_seq_len);
    for case in &GLOBAL_CASES {
        let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
        let scoring = SequenceSimilarityScore::for_strings(&case.query, &case.template, blosum62);
        let score = aligner.align(&scoring, -10, -2);
        let path = aligner.backtrace();
        let (ali_q, ali_t) = aligned_strings(&path, &case.query, &case.template, '-');
        assert_eq!(score, case.score);
        assert_eq!(aligner.recent_score(), case.score);
        assert_eq!(ali_q, case.aligned_query);
        assert_eq!(ali_t, case.aligned_template);
        assert_eq!(path.to_string(), case.alignment);

        let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
        let scoring = SequenceSimilarityScore::for_strings(&case.template,&case.query,  blosum62);
        let score = aligner.align(&scoring, -10, -2);
        let path = aligner.backtrace();
        let (ali_q, ali_t) = aligned_strings(&path, &case.template, &case.query, '-');
        assert_eq!(score, case.score);
        assert_eq!(aligner.recent_score(), case.score);
        assert_eq!(ali_t, case.aligned_query);
        assert_eq!(ali_q, case.aligned_template);
    }
}
