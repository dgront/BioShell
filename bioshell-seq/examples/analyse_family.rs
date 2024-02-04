use std::collections::HashMap;
use std::env;
use clap::{Parser};
use log::{info};

use bioshell_seq::sequence::{Sequence, FastaIterator};
use bioshell_io::open_file;
use bioshell_seq::alignment::{GlobalAligner, aligned_sequences, AlignmentReporter, AlignmentStatistics};
use bioshell_seq::scoring::{SequenceSimilarityScore, SubstitutionMatrixList};
use bioshell_statistics::Histogram;

#[derive(Parser, Debug)]
#[clap(name = "analyse_family")]
#[clap(about = "Calculates global sequence alignment of amino acid sequences", long_about = None)]
struct Args {
    /// query sequence(s): either a FASTA string or a name of a file in FASTA format
    #[clap(long, short='q')]
    query: String,
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
        let reader = open_file(&seq_or_fname);
        let seq_iter = FastaIterator::new(reader);
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
            if lower_triangle && template==query { break }
            scoring.query_from_sequence(query);
            aligner.align(&scoring, gap_open, gap_extend);
            let path = aligner.backtrace();
            let (ali_q, ali_t) = aligned_sequences(&path, &query, &template, '-');
            reporter.report(&ali_q, &ali_t);
        }
    }
}

struct SimilarityHistogramByQuery {
    histogram_bin_width: f64,
    similarities: HashMap<String, Histogram>
}

impl SimilarityHistogramByQuery {
    pub fn new(histogram_bin_width: f64) -> SimilarityHistogramByQuery {
        SimilarityHistogramByQuery { histogram_bin_width, similarities: Default::default() }
    }
}

impl AlignmentReporter for SimilarityHistogramByQuery {
    fn report(&mut self, aligned_query: &Sequence, aligned_template: &Sequence) {
        let stats = AlignmentStatistics::from_sequences(aligned_query, aligned_template, 32);

        if !self.similarities.contains_key(&stats.query_header) {
            self.similarities.insert(stats.query_header.clone(), Histogram::by_bin_width(self.histogram_bin_width));
        } else {
            let h = self.similarities.get_mut(&stats.query_header).unwrap();
            h.insert(stats.percent_identity());
        }
    }
}

impl Drop for SimilarityHistogramByQuery {
    fn drop(&mut self) {
        for (key, histogram) in &self.similarities {
            println!("{} min: {}, max: {}, mode: {}\n{}", key,
                     histogram.min().unwrap() as f64 * self.histogram_bin_width,
                     histogram.max().unwrap() as f64 * self.histogram_bin_width,
                     histogram.tallest().unwrap() as f64 * self.histogram_bin_width,
                     histogram.draw_horizonaly(0.0, 100.0, 10));
        }
    }
}


pub fn main() {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();

    let queries = get_sequences(&args.query, "query");

    align_all_pairs(&queries, &queries, SubstitutionMatrixList::BLOSUM62,
                    args.open, args.extend, true, &mut SimilarityHistogramByQuery::new(2.5));
}