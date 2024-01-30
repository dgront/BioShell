use std::fmt::{Display, Formatter};
use crate::sequence::{count_identical, len_ungapped, Sequence};

/// Provides staple statistics for a sequence alignment
pub struct  AlignmentStatistics {
    /// identifies the query sequence
    pub query_header: String,
    /// identifies the template sequence
    pub template_header: String,
    /// counts identical positions within the alignment
    pub n_identical: usize,
    /// length of the query sequence without gaps
    pub query_length: usize,
    /// length of the template sequence without gaps
    pub template_length: usize,
}

impl AlignmentStatistics {
    /// Creates the [AlignmentStatistics] for a given pair of aligned sequences
    pub fn from_sequences(aligned_query: &Sequence, aligned_template: &Sequence, header_length: usize) -> AlignmentStatistics {
        let query_header = aligned_query.description()[0..header_length].to_string();
        let template_header = aligned_template.description()[0..header_length].to_string();
        let n_identical = count_identical(aligned_query, aligned_template).unwrap();
        let query_length = len_ungapped(aligned_query);
        let template_length = len_ungapped(aligned_template);

        AlignmentStatistics {
            query_header, template_header,
            n_identical, query_length, template_length,
        }
    }

    /// Computes the sequence identity
    pub fn percent_identity(&self) -> f64 {
        self.n_identical as f64 / self.query_length.min(self.template_length) as f64 * 100.0
    }
}

impl Display for AlignmentStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:35} {:35} {:6.2} % {:3} {:4} {:4}", self.query_header, self.template_header,
               self.percent_identity(), self.n_identical, self.query_length, self.template_length)
    }
}