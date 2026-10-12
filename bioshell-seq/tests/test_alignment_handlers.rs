#[cfg(test)]
mod tests {
    use bioshell_seq::alignment::{format_pairwise, LocalAlignment};
    use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
    use bioshell_seq::sequence::Sequence;

    #[test]
    fn test_pairwise_formatting() {
        let query = Sequence::new("query", "MNPQRAAVVGSTKLTQAAIDYKGYEVNVRHDF");
        let template = Sequence::new("template", "WCEGAVVGSTKLSQAAIKGYEVTVKLPF");
        let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
        let scoring = SequenceSimilarityScore::for_sequences(&query, &template, blosum62);
        let mut aligner = LocalAlignment::new(33);
        let _score = aligner.align(&scoring, -10, -2);
        let path = aligner.backtrace();

        let out = format_pairwise(&query, &template, &path, 10, 20);
        let expected = "query           7 AVVGSTKLTQAAIDYKGYEV 26
                  |||||||| ||||  |||||
template        5 AVVGSTKLSQAAI--KGYEV 22

query          27 NVRHDF 32
                   |   |
template       23 TVKLPF 28";
        assert_eq!(out.trim(), expected.trim());
    }
}