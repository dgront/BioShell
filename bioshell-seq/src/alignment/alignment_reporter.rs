use crate::alignment::AlignmentStatistics;
use crate::sequence::{len_ungapped_str, Sequence};

/// Reports a sequence alignment calculated by a sequence alignment algorithm.
pub trait AlignmentReporter {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence);
}

/// Reporter that combines multiple reporters.
pub struct MultiReporter {
    reporters: Vec<Box<dyn AlignmentReporter>>,
}
impl MultiReporter {
    pub fn new() -> Self { MultiReporter { reporters: vec![] } }
    pub fn add_reporter(&mut self, reporter: Box<dyn AlignmentReporter>) {
        self.reporters.push(reporter);
    }
    pub fn count_reporters(&self) -> usize { self.reporters.len() }
}

impl AlignmentReporter for MultiReporter {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        for reporter in &mut self.reporters {
            reporter.report(aligned_query, aligned_template);
        }
    }
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
/// let template = Sequence::from_str("template", "ALRIV");
/// let mut reporter = PrintAsPairwise::new(5, 10);
/// reporter.report(&query, &template);
/// ```
pub struct PrintAsPairwise {
    /// maximum number of characters available to print a sequence name
    pub seq_name_width: usize,
    pub alignment_width: usize
}

impl PrintAsPairwise {
    pub fn new(seq_name_width: usize, alignment_width: usize) -> PrintAsPairwise {
        PrintAsPairwise { seq_name_width, alignment_width }
    }
}

impl AlignmentReporter for PrintAsPairwise {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let q_name = format!("{:len$}", aligned_query.description_n(self.seq_name_width), len = self.seq_name_width);
        let t_name = format!("{:len$}", aligned_template.description_n(self.seq_name_width), len = self.seq_name_width);
        let query_chunks: Vec<String> = aligned_query.to_string(0).chars().collect::<Vec<_>>()
                .chunks(self.alignment_width).map(|chunk| chunk.into_iter().collect()).collect();
        let template_chunks: Vec<String> = aligned_template.to_string(0).chars().collect::<Vec<_>>()
            .chunks(self.alignment_width).map(|chunk| chunk.into_iter().collect()).collect();

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