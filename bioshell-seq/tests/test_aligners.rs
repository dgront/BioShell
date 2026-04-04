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

use bioshell_seq::alignment::{LocalAlignment};

struct LocalAlignmentTestCase {
    query: &'static str,
    template: &'static str,
    aligned_query: &'static str,
    aligned_template: &'static str,
    alignment: &'static str,
    query_start: usize,
    template_start: usize,
    score: i32,
}

static LOCAL_CASES: [LocalAlignmentTestCase; 3] = [
    LocalAlignmentTestCase {
        query: "MNPQRAAVVGSTKLTQAAIDYKGYEVNVRHDF",
        template: "WCEGAVVGSTKLSQAAIKGYEVTVKLPF",
        aligned_query: "AVVGSTKLTQAAIDYKGYEVNVRHDF",
        aligned_template: "AVVGSTKLSQAAI--KGYEVTVKLPF",
        alignment:        "*************||***********",
        query_start: 6,
        template_start: 4,
        score: 77,
    },
    LocalAlignmentTestCase {
        query: "GHMTRQILDSLGNPTAEEVKAAFDKLMNQERWP",
        template: "PPYQILDSLGAEEVRAAFDKLQERFA",
        aligned_query: "QILDSLGNPTAEEVKAAFDKLMNQERW",
        aligned_template: "QILDSLG---AEEVRAAFDKL--QERF",
        alignment:        "*******|||***********||****",
        query_start: 5,
        template_start: 3,
        score: 72,
    },
    LocalAlignmentTestCase {
        query: "MAVR",
        template: "QQMAVRPP",
        aligned_query: "MAVR",
        aligned_template: "MAVR",
        alignment: "****",
        query_start: 0,
        template_start: 2,
        score: 18,
    },
];

#[test]
fn test_local_alignment() {
    let max_seq_len = LOCAL_CASES.iter().map(|c| c.query.len().max(c.template.len())).max().unwrap();
    let mut aligner = LocalAlignment::new(max_seq_len);

    for case in &LOCAL_CASES {
        let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
        let scoring = SequenceSimilarityScore::for_strings(case.query, case.template, blosum62);

        let score = aligner.align(&scoring, -10, -2);
        let (path, query_start, template_start) = aligner.backtrace();
        let (ali_q, ali_t) = aligned_strings(
            &path,
            &case.query[query_start..],
            &case.template[template_start..],
            '-',
        );

        assert_eq!(score, case.score);
        assert_eq!(aligner.recent_score(), case.score);
        assert_eq!(query_start, case.query_start);
        assert_eq!(template_start, case.template_start);
        assert_eq!(ali_q, case.aligned_query);
        assert_eq!(ali_t, case.aligned_template);
        assert_eq!(path.to_string(), case.alignment);

        let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
        let scoring = SequenceSimilarityScore::for_strings(case.template, case.query, blosum62);

        let score = aligner.align(&scoring, -10, -2);
        let (path, query_start, template_start) = aligner.backtrace();
        let (ali_q, ali_t) = aligned_strings(
            &path,
            &case.template[query_start..],
            &case.query[template_start..],
            '-',
        );

        assert_eq!(score, case.score);
        assert_eq!(aligner.recent_score(), case.score);
        assert_eq!(ali_q, case.aligned_template);
        assert_eq!(ali_t, case.aligned_query);
    }
}
