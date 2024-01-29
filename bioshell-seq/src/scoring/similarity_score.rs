use crate::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
use crate::sequence::Sequence;

pub trait SimilarityScore {
    fn score(&self, i_pos: usize, j_pos: usize) -> i32;
    fn template_length(&self) -> usize;
    fn query_length(&self) -> usize;
}

pub struct SequenceSimilarityScore {
    query: Vec<u8>,
    template: Vec<u8>,
    matrix: SubstitutionMatrix
}

impl SequenceSimilarityScore {
    pub fn new(matrix: SubstitutionMatrixList) -> SequenceSimilarityScore {
        SequenceSimilarityScore{
            query: vec![], template: vec![], matrix:SubstitutionMatrix::load(matrix),
        }
    }
    pub fn for_sequences(query: &Sequence, template: &Sequence, matrix: SubstitutionMatrix) -> SequenceSimilarityScore {
        SequenceSimilarityScore{
            query: query.as_u8().iter().map(|a|matrix.aa_index(*a)).collect(),
            template: template.as_u8().iter().map(|a|matrix.aa_index(*a)).collect(),
            matrix,
        }
    }

    pub fn for_strings(query: &str, template: &str, matrix: SubstitutionMatrix) -> SequenceSimilarityScore {
        let query = Sequence::from_str("query", query);
        let template = Sequence::from_str("template", template);

        SequenceSimilarityScore::for_sequences(&query, &template, matrix)
    }

    pub fn query_from_str(&mut self, query: &str) {
        self.query = query.as_bytes().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    pub fn template_from_str(&mut self, template: &str) {
        self.template = template.as_bytes().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    pub fn query_from_sequence(&mut self, query: &Sequence) {
        self.query = query.as_u8().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    pub fn template_from_sequence(&mut self, template: &Sequence) {
        self.template = template.as_u8().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }
}

impl SimilarityScore for SequenceSimilarityScore {
    fn score(&self, q_pos: usize, t_pos: usize) -> i32 {
        self.matrix.score_by_index(self.query[q_pos], self.template[t_pos])
    }

    fn template_length(&self) -> usize { self.template.len() }

    fn query_length(&self) -> usize { self.query.len() }

}