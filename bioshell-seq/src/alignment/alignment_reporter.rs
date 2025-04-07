use std::collections::HashMap;
use std::fs::File;
use bioshell_io::out_writer;
use crate::alignment::AlignmentStatistics;
use crate::sequence::{parse_sequence_id, len_ungapped, len_ungapped_str, Sequence};

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
    pub header_width: usize,
    pub infer_seq_id: bool
}

impl SimilarityReport {
    pub fn new(header_width: usize, infer_seq_id: bool) -> SimilarityReport {
        SimilarityReport{ header_width, infer_seq_id }
    }
}

impl Default for SimilarityReport {
    fn default() -> Self { SimilarityReport::new(32, false) }
}

impl AlignmentReporter for SimilarityReport {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let mut stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, self.header_width);
        if self.infer_seq_id {
            stats.query_header = parse_sequence_id(&stats.query_header).to_string();
            stats.template_header = parse_sequence_id(&stats.template_header).to_string();
        }
        println!("{}", stats);
    }
}

/// Prints sequence identity matrix for a number of sequences aligned with one another.
///
/// # Example
/// ```
/// use bioshell_seq::alignment::{AlignmentReporter, IdentityMatrixReporter};
/// use bioshell_seq::sequence::Sequence;
/// let query = Sequence::from_str("query", "AL-IV");
/// let template = Sequence::from_str("tmplt", "ALRIV");
/// let mut reporter = IdentityMatrixReporter::new(5, false, "stdout");
/// reporter.report(&query, &template);
/// assert_eq!(reporter.sequence_index("query"), Some(0));
/// assert_eq!(reporter.sequence_index("tmplt"), Some(1));
/// assert_eq!(reporter.n_identical_residues(0, 1), 4);
pub struct IdentityMatrixReporter {
    pub header_width: usize,
    pub out_fname: String,
    pub infer_seq_id: bool,
    identity_matrix: Vec<Vec<usize>>,
    sequence_order: HashMap<String, usize>,
}

impl IdentityMatrixReporter {
    pub fn new(header_width: usize, infer_seq_id: bool, out_fname: &str) -> IdentityMatrixReporter {
        IdentityMatrixReporter {
            header_width, out_fname: out_fname.to_string(), infer_seq_id,
            identity_matrix: vec![], sequence_order: HashMap::new()
        }
    }

    /// Returns the index of the sequence in the identity matrix.
    ///
    /// Sequences have been recorded in the order they were added to this alignment reporter;
    /// the query sequence is inserted before the template.
    pub fn sequence_index(&mut self, seq_name: &str) -> Option<usize> {
        if !self.sequence_order.contains_key(seq_name) {
           return None;
        }
        Some(self.sequence_order[seq_name])
    }

    /// Returns the number of sequences reported so far
    pub fn num_sequences(&self) -> usize { self.sequence_order.len() }

    pub fn n_identical_residues(&self, query_idx: usize, tmplt_idx: usize) -> usize {
        if query_idx < tmplt_idx {
            self.identity_matrix[tmplt_idx][query_idx]
        } else {
            self.identity_matrix[query_idx][tmplt_idx]
        }
    }
}

impl Drop for IdentityMatrixReporter {
    fn drop(&mut self) {
        // ---------- Open the file for writing
        let mut file = out_writer(&self.out_fname, false);
        let err_msg = format!("Failed to write to file: {}", self.out_fname);

        // ---------- Sequences in the order they were added to the alignment reporter
        let mut sorted_keys: Vec<(&String, &usize)> = self.sequence_order.iter().collect();
        sorted_keys.sort_by_key(|&(_, &index)| index);

        // ---------- Write all sequence headers
        for (key, _) in &sorted_keys {
            writeln!(file, "{}", key).expect(&err_msg);
        }
        writeln!(file).expect(&err_msg); // Add a blank line to separate blocks

        // Write the header and corresponding matrix values (second block)
        for (key, &idx) in &sorted_keys {
            // if requested, extract the seq-id from a sequence header
            let mut truncated_key = if self.infer_seq_id { parse_sequence_id(key).to_string()} else { (*key).clone() };
            // Write the truncated key (header_width characters)
            truncated_key = truncated_key.split_whitespace().next().unwrap().chars().take(self.header_width).collect();
            write!(file, "{:<width$}", truncated_key, width = self.header_width).expect(&err_msg);

            // Write the values from identity_matrix[idx]
            if let Some(row) = self.identity_matrix.get(idx) {
                for &value in row { write!(file, " {:3}", value).expect(&err_msg); }
            }

            writeln!(file).expect(&err_msg);
        }
    }
}

impl AlignmentReporter for IdentityMatrixReporter {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let q_name = aligned_query.description();
        let t_name = aligned_template.description();
        if !self.sequence_order.contains_key(q_name) {
            self.sequence_order.insert(q_name.to_string(), self.identity_matrix.len());
            self.identity_matrix.push(vec![0; self.sequence_order.len()]);
            let idx = self.sequence_order.len() - 1;
            self.identity_matrix[idx][idx] = len_ungapped(aligned_query);
        }
        if !self.sequence_order.contains_key(t_name) {
            self.sequence_order.insert(t_name.to_string(), self.identity_matrix.len());
            self.identity_matrix.push(vec![0; self.sequence_order.len()]);
            let idx = self.sequence_order.len() - 1;
            self.identity_matrix[idx][idx] = len_ungapped(aligned_template);
        }
        let q_idx = self.sequence_order[q_name];
        let t_idx = self.sequence_order[t_name];
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, self.header_width);
        if t_idx > q_idx {
            self.identity_matrix[t_idx][q_idx] = stats.n_identical;
        } else {
            self.identity_matrix[q_idx][t_idx] = stats.n_identical;
        }
    }
}