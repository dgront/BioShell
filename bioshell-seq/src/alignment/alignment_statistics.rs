use std::fmt::{Display, Formatter};
use crate::sequence::{count_identical, len_ungapped, Sequence};

/// Provides staple statistics for a sequence alignment.
///
/// The struct is created based on two [`Sequence`]s that are assumed to be aligned. It provides
/// description line for both sequences (trimmed to desired length) and numerical statistics computed from
/// the alignment such as the number of identical residues, the length of the alignment, etc.
///
/// # Example
/// ```
/// use bioshell_seq::alignment::AlignmentStatistics;
/// let aligned_query = "EIIIDSYNQFSDR----SYQFMTPSLFVR";
/// let aligned_tmplt = "ETVKEAYDLYPDRRYFGSFQFLYPSLFLR";
/// let stats = AlignmentStatistics::from_strings("query", aligned_query, "template", aligned_tmplt, 5);
/// let stats_printed = format!("{}", stats);
/// assert_eq!(stats.query_length, 25);
/// assert_eq!(stats.template_length, 29);
/// assert_eq!(stats.n_identical, 12);
/// assert_eq!(stats_printed, "query templ  48.00 %  12   25   29".to_string());
/// ```
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
    /// maximum number of characters a sequence name can take
    header_length: usize
}

impl AlignmentStatistics {
    pub fn from_strings(query_name: &str, query_sequence: &str,
        template_name: &str, template_sequence: &str, name_width: usize) -> AlignmentStatistics {
        let q = Sequence::from_str(query_name, query_sequence);
        let t = Sequence::from_str(template_name, template_sequence);

        return AlignmentStatistics::from_sequences(&q, &t, name_width);
    }
    /// Creates the [AlignmentStatistics] for a given pair of aligned sequences

    pub fn from_sequences(aligned_query: &Sequence, aligned_template: &Sequence, header_length: usize) -> AlignmentStatistics {
        let query_header = aligned_query.description_n(header_length);
        let template_header = aligned_template.description_n(header_length);
        let n_identical = count_identical(aligned_query, aligned_template).unwrap();
        let query_length = len_ungapped(aligned_query);
        let template_length = len_ungapped(aligned_template);

        AlignmentStatistics {
            query_header: query_header.to_string(), template_header: template_header.to_string(),
            n_identical, query_length, template_length, header_length
        }
    }

    /// Computes the sequence identity
    pub fn percent_identity(&self) -> f64 {
        self.n_identical as f64 / self.query_length.min(self.template_length) as f64 * 100.0
    }
}

impl Display for AlignmentStatistics {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:len$} {:len$} {:6.2} % {:3} {:4} {:4}", self.query_header, self.template_header,
               self.percent_identity(), self.n_identical, self.query_length, self.template_length,
               len = self.header_length)
    }
}