use crate::scoring::{SubstitutionMatrix, SubstitutionMatrixList};
use crate::sequence::Sequence;

/// Trait for scoring the similarity of two sequences.
///
/// This is used by the sequence alignment algorithms to calculate the score of a particular alignment.
/// The score is calculated by the `score(i, j)` method, which takes the positions of the two sequences
/// and returns the respective score.
pub trait SimilarityScore {
    /// Calculates the score of the similarity of two sequences at the given positions.
    fn score(&self, i_pos: usize, j_pos: usize) -> i32;
    /// Returns the length of the template sequence.
    fn template_length(&self) -> usize;
    /// Returns the length of the query sequence.
    fn query_length(&self) -> usize;
    /// Returns `true` if the characters at the given positions are identical, `false` otherwise.
    fn is_identity(&self, _i: usize, _j: usize) -> bool;
}

/// Calculates the [`SimilarityScore`] of two sequences using a [`SubstitutionMatrix`].
///
/// The query and template sequences are stored as vectors of indices corresponding to the characters
/// in the substitution matrix. The score is calculated with [`score()`](score) by looking up the score for the respective indices
/// in the substitution matrix. Either the query or the template sequence can be replaced by a new one
/// using the [`query_from_str()`](query_from_str), [`template_from_str()`](template_from_str),
/// [`query_from_sequence()`](query_from_sequence) and [`template_from_sequence()`](template_from_sequence) methods.
///
/// # Example
///```
/// use bioshell_seq::scoring::{SequenceSimilarityScore, SimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
/// use bioshell_seq::sequence::Sequence;
/// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
/// let mut similarity_score = SequenceSimilarityScore::for_strings("AKLIY", "AGKIY", blosum62);
/// assert_eq!(similarity_score.score(1, 1), -2); // K-G is -2 by BLOSUM62
/// similarity_score.template_from_str("GKLLY");
/// assert_eq!(similarity_score.score(1, 1), 5); // K-K is +5 by BLOSUM62
/// ```
pub struct SequenceSimilarityScore {
    query: Vec<u8>,
    template: Vec<u8>,
    matrix: SubstitutionMatrix
}

impl SequenceSimilarityScore {
    /// Creates a new [`SequenceSimilarityScore`] instance with the specified substitution matrix.
    ///
    /// The query and template sequences are initialized as empty vectors and must be replaced later.
    /// This constructor is useful when the same substitution matrix is used to calculate the similarity score for multiple pairs of sequences.
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::scoring::{SequenceSimilarityScore, SimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
    /// let sequences = vec!["AKLIY", "AGKIY", "GKLLY"];
    /// # let mut observed_values: Vec<i32> = Vec::new();
    /// let mut similarity = SequenceSimilarityScore::new(SubstitutionMatrixList::BLOSUM62);
    /// for i in 0..sequences.len() {
    ///     similarity.query_from_str(sequences[i]);
    ///     for j in (i + 1)..sequences.len() {
    ///         similarity.template_from_str(sequences[j]);
    ///         let sim: i32 = (0..5).map(|i| similarity.score(i, i)).sum();
    /// #       observed_values.push(sim);
    ///         println!("Similarity of seq {} and seq {} is {}", i, j, sim);
    ///     }
    /// }
    /// # assert_eq!(observed_values, vec![11, 18, 5]);
    /// ```
    pub fn new(matrix: SubstitutionMatrixList) -> SequenceSimilarityScore {
        SequenceSimilarityScore{
            query: vec![], template: vec![], matrix:SubstitutionMatrix::load(matrix),
        }
    }

    /// Creates a new [`SequenceSimilarityScore`] instance for the given query and template sequences and substitution matrix.
    ///
    /// ```
    /// use bioshell_seq::scoring::{SequenceSimilarityScore, SimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
    /// use bioshell_seq::sequence::Sequence;
    /// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
    /// let query = Sequence::from_str("query", "AKLIY");
    /// let template = Sequence::from_str("template", "AGKIY");
    /// let similarity_score = SequenceSimilarityScore::for_sequences(&query, &template, blosum62);
    /// assert_eq!(similarity_score.score(1, 1), -2); // K-G is -2 by BLOSUM62
    /// ```
    pub fn for_sequences(query: &Sequence, template: &Sequence, matrix: SubstitutionMatrix) -> SequenceSimilarityScore {
        SequenceSimilarityScore{
            query: query.as_u8().iter().map(|a|matrix.aa_index(*a)).collect(),
            template: template.as_u8().iter().map(|a|matrix.aa_index(*a)).collect(),
            matrix,
        }
    }

    /// Creates a new [`SequenceSimilarityScore`] instance for the given query and template sequence strings and substitution matrix.
    ///
    /// ```
    /// use bioshell_seq::scoring::{SequenceSimilarityScore, SimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};
    /// use bioshell_seq::sequence::Sequence;
    /// let blosum62 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM62);
    /// let similarity_score = SequenceSimilarityScore::for_strings("AKLIY", "AGKIY", blosum62);
    /// assert_eq!(similarity_score.score(1, 1), -2); // K-G is -2 by BLOSUM62
    /// ```
    pub fn for_strings(query: &str, template: &str, matrix: SubstitutionMatrix) -> SequenceSimilarityScore {
        let query = Sequence::from_str("query", query);
        let template = Sequence::from_str("template", template);

        SequenceSimilarityScore::for_sequences(&query, &template, matrix)
    }

    /// Replaces the query sequence with a new one given as a string.
    pub fn query_from_str(&mut self, query: &str) {
        self.query = query.as_bytes().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    /// Replaces the template sequence with a new one given as a string.
    pub fn template_from_str(&mut self, template: &str) {
        self.template = template.as_bytes().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    /// Replaces the query sequence with a new one given as a [`Sequence`].
    pub fn query_from_sequence(&mut self, query: &Sequence) {
        self.query = query.as_u8().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }

    /// Replaces the template sequence with a new one given as a [`Sequence`].
    pub fn template_from_sequence(&mut self, template: &Sequence) {
        self.template = template.as_u8().iter().map(|a| self.matrix.aa_index(*a)).collect();
    }
}

impl SimilarityScore for SequenceSimilarityScore {
    #[inline(always)]
    fn score(&self, q_pos: usize, t_pos: usize) -> i32 {
        self.matrix.score_by_index(self.query[q_pos], self.template[t_pos])
    }

    #[inline(always)]
    fn template_length(&self) -> usize { self.template.len() }

    #[inline(always)]
    fn query_length(&self) -> usize { self.query.len() }

    #[inline(always)]
    fn is_identity(&self, i: usize, j: usize) -> bool {
        self.query[i] == self.template[j]
    }
}