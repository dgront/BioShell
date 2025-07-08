use std::collections::HashMap;
use std::env;
use std::io::Write;
use clap::{Parser};
#[allow(unused_imports)]
use log::{info};

use bioshell_seq::sequence::{Sequence, load_sequences};
use bioshell_io::{out_writer};
use bioshell_seq::alignment::{AlignmentReporter, AlignmentStatistics, align_all_pairs, PrintAsPairwise, SimilarityReport, MultiReporter};
use bioshell_seq::scoring::{SubstitutionMatrixList};
use bioshell_seq::SequenceError;
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
    extend: i32,
    /// print pairwise alignments for every pair of aligned sequences
    #[clap(long, action)]
    pairwise: bool,
    /// print a histogram of sequence identity values for each sequence
    #[clap(long)]
    histograms: Option<String>,
    /// print sequence identity report (default)
    #[clap(long, action)]
    report: bool,
    /// length of a sequence name to print; longer names will be cut that size
    #[clap(long, short='w', default_value = "20")]
    name_width: usize,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

struct SimilarityHistogramByQuery {
    histogram_bin_width: f64,
    output_file_name: String,
    similarities: HashMap<String, Histogram>
}

impl SimilarityHistogramByQuery {
    pub fn new(histogram_bin_width: f64, output_file_name: &str) -> SimilarityHistogramByQuery {
        SimilarityHistogramByQuery { histogram_bin_width, output_file_name: output_file_name.to_string(),
            similarities: Default::default()
        }
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
        let mut writer = out_writer(&self.output_file_name, false);
        for (key, histogram) in &self.similarities {
            writer.write(format!("{} min: {}, max: {}, mode: {}\n{}", key,
                                 histogram.min().unwrap() as f64 * self.histogram_bin_width,
                                 histogram.max().unwrap() as f64 * self.histogram_bin_width,
                                 histogram.tallest().unwrap() as f64 * self.histogram_bin_width,
                                 histogram.draw_horizonaly(0.0, 100.0, 10)).as_bytes()).unwrap();
        }
    }
}


pub fn main() -> Result<(), SequenceError> {
    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    let build_time = env!("BUILD_TIME");
    let git_commit_md5 = env!("GIT_COMMIT_MD5");

    info!("Build time: {}", build_time);
    info!("Git commit MD5 sum: {}", git_commit_md5);

    let name_width = args.name_width;

    let queries = load_sequences(&args.query, "query")?;

    let mut multireports = MultiReporter::new();
    if let Some(fname) = args.histograms {
        multireports.add_reporter(Box::new(SimilarityHistogramByQuery::new(2.5, &fname)));
    }
    if args.pairwise { multireports.add_reporter(Box::new(PrintAsPairwise::new(name_width, 80))); }
    if args.report { multireports.add_reporter(Box::new(SimilarityReport::new(name_width, false))); }
    if multireports.count_reporters() == 0 {  multireports.add_reporter(Box::new(SimilarityReport::new(name_width, false))); }

    align_all_pairs(&queries, &queries, SubstitutionMatrixList::BLOSUM62,
                    args.open, args.extend, true, &mut multireports);

    Ok(())
}