use std::collections::HashMap;
use std::env;
use std::fmt::{Display, Formatter};
use clap::{Parser};
use log::{info};

use bioshell_seq::sequence::{Sequence, FastaIterator, count_identical, len_ungapped};
use bioshell_io::open_file;
use bioshell_seq::alignment::{GlobalAligner, aligned_sequences};
use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrix, SubstitutionMatrixList};

#[derive(Parser, Debug)]
#[clap(name = "needleman_wunsh")]
#[clap(about = "Calculates global sequence alignment of amino acid sequences", long_about = None)]
struct Args {
    /// query sequence(s): either a FASTA string or a name of a file in FASTA format
    #[clap(long, short='q')]
    query: String,
    /// template sequence(s): either a FASTA string or a name of a file in FASTA format
    #[clap(long, short='t')]
    template: Option<String>,
    /// gap opening penalty
    #[clap(long, default_value = "-10", short='o')]
    open: i32,
    /// gap extension penalty
    #[clap(long, default_value = "-2", short='e')]
    extend: i32
}

/// Returns a list of Sequences for a given input string.
///
/// If the input string contains a dot character, it's assumed to be a file name. This function
/// attempts to open that file as FASTA and load all sequences stored.
/// Otherwise it's assumed the string is an amino acid sequence by itself; it's converted into sequence and returned
/// as the only element of a vector.
fn get_sequences(seq_or_fname: &String, seq_name: &str) -> Vec<Sequence> {
    if seq_or_fname.contains(".") {
        let mut out : Vec<Sequence> = vec![];
        let reader = open_file(&seq_or_fname);
        let mut seq_iter = FastaIterator::new(reader);
        return  seq_iter.collect();
    }

    return vec![Sequence::from_str(seq_name, seq_or_fname)];
}

fn align_all_pairs<R: AlignmentReporter>(queries: &Vec<Sequence>, templates: &Vec<Sequence>,
        matrix: SubstitutionMatrixList, gap_open: i32, gap_extend: i32, lower_triangle: bool,
        reporter: &mut R) {

    let mut max_length = queries.iter().map(|s| s.len()).max().unwrap();
    max_length = max_length.max(templates.iter().map(|s| s.len()).max().unwrap());
    let mut aligner = GlobalAligner::new(max_length);
    let mut scoring = SequenceSimilarityScore::new( matrix);

    for template in templates {
        scoring.template_from_sequence(template);
        for query in queries {
            if lower_triangle && template==query { continue }
            scoring.query_from_sequence(query);
            aligner.align(&scoring, gap_open, gap_extend);
            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, &query, &template, '-');
            reporter.report(&ali_q, &ali_t);
        }
    }
}

struct  AlignmentStatistics {
    pub query_header: String,
    pub template_header: String,
    pub n_identical: usize,
    pub query_length: usize,
    pub template_length: usize,
}

impl AlignmentStatistics {
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

trait AlignmentReporter {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence);
}

struct PrintAsFasta;

impl AlignmentReporter for PrintAsFasta {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        println!("{}\n{}", aligned_query, aligned_template);
    }
}

struct SimilarityReport;

impl AlignmentReporter for SimilarityReport {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, 32);
        println!("{}", stats);
    }
}


struct SimilarityByQuery {
    key_length: usize,
    min_max: HashMap<String, (f64, f64)>
}

impl SimilarityByQuery {
    pub fn new(key_length: usize) -> SimilarityByQuery {
        SimilarityByQuery{ key_length: key_length, min_max: Default::default() }
    }
}

impl AlignmentReporter for SimilarityByQuery {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, 32);

        if !self.min_max.contains_key(&stats.query_header) {
            self.min_max.insert(stats.query_header.clone(), (100.0, 0.0));
        } else {
            let vals = self.min_max.get_mut(&stats.query_header).unwrap();
            vals.0 = vals.0.min(stats.percent_identity());
            vals.1 = vals.1.max(stats.percent_identity());
        }
    }
}

impl Drop for SimilarityByQuery {
    fn drop(&mut self) {
        for (key, val) in &self.min_max {
            println!("{} {:6.2} {:6.2}", key, val.0, val.1);
        }
    }
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();

    let queries = get_sequences(&args.query, "query");

    if let Some(tmpl) = args.template {
        let templates = get_sequences(&tmpl, "template");
        align_all_pairs(&queries, &templates, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, true, &mut SimilarityReport);

    } else {
        align_all_pairs(&queries, &queries, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, false, &mut SimilarityByQuery::new(32));
    }

}