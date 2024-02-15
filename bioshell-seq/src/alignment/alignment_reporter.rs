use crate::alignment::AlignmentStatistics;
use crate::sequence::{len_ungapped_str, Sequence};

/// Reports a sequence alignment
pub trait AlignmentReporter {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence);
}

/// Prints a sequence alignment in the FASTA format
pub struct PrintAsFasta;

impl AlignmentReporter for PrintAsFasta {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        println!("{}\n{}", aligned_query, aligned_template);
    }
}

/// Prints a sequence alignment in the *pairwise* format
///
/// # Example
/// ```
/// use bioshell_seq::alignment::{AlignmentReporter, PrintAsPairwise};
/// use bioshell_seq::sequence::Sequence;
/// let query = Sequence::from_str("query", "AL-IV");
/// let template = Sequence::from_str("query", "ALRIV");
/// let mut reporter = PrintAsPairwise::new(10);
/// reporter.report(&query, &template);
/// ```
pub struct PrintAsPairwise {
    /// maximum number of characters available to print a sequence name
    ///
    /// By default this value is set 16
    pub seq_name_width: usize,
    width: usize
}

impl PrintAsPairwise {
    pub fn new(output_width: usize) -> PrintAsPairwise { PrintAsPairwise { seq_name_width: 16, width: output_width } }
}


impl AlignmentReporter for PrintAsPairwise {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let len = aligned_query.description().len().max(aligned_template.description().len()).min(self.seq_name_width);
        let q_name = format!("{:len$}", aligned_query.description(), len = len);
        let t_name = format!("{:len$}", aligned_template.description(), len = len);
        let query_chunks: Vec<String> = aligned_query.to_string().chars().collect::<Vec<_>>()
                .chunks(self.width).map(|chunk| chunk.into_iter().collect()).collect();
        let template_chunks: Vec<String> = aligned_template.to_string().chars().collect::<Vec<_>>()
            .chunks(self.width).map(|chunk| chunk.into_iter().collect()).collect();

        let mut q_from = 1;
        let mut t_from = 1;
        let numeric_width = 5;
        let num_spacer = std::iter::repeat(" ").take(numeric_width).collect::<String>();
        for (q, t) in query_chunks.iter().zip(&template_chunks) {
            let iter1 = q.chars();
            let iter2 = t.chars();
            let mid_name = std::iter::repeat(" ").take(q_name.len()).collect::<String>();
            let middle: Vec<u8> = iter1.zip(iter2).map(|(c1, c2)| if c1==c2 { b'|'} else {b' '} ).collect();
            let q_add = len_ungapped_str(q);
            let t_add = len_ungapped_str(t);
            println!("{} {:width$} {} {:width$}\n{} {} {}\n{} {:width$} {} {:width$}", q_name, q_from, q, q_from + q_add,
                     mid_name, num_spacer, String::from_utf8(middle).unwrap(),
                     t_name, t_from, t, t_from + t_add, width = numeric_width);
            q_from += q_add - 1;
            t_from += t_add - 1;
        }
    }
}

/// Prints staple statistics for a given alignment
pub struct SimilarityReport {
    pub header_width: usize
}

impl SimilarityReport {
    pub fn new(header_width: usize) -> SimilarityReport { SimilarityReport{ header_width } }
}

impl Default for SimilarityReport {
    fn default() -> Self { SimilarityReport::new(32) }
}

impl AlignmentReporter for SimilarityReport {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, self.header_width);
        println!("{}", stats);
    }
}