use std::env;
use std::fmt::{Display};
use clap::{Parser};
use log::{info};

use bioshell_seq::sequence::{Sequence, FastaIterator};
use bioshell_io::open_file;
use bioshell_seq::alignment::{PrintAsPairwise, AlignmentReporter, SimilarityReport, align_all_pairs};
use bioshell_seq::scoring::{SubstitutionMatrixList};

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
    extend: i32,
    /// print pairwise alignments for every pair of aligned sequences
    #[clap(long, action)]
    pairwise: bool,
    /// print sequence identity report (default)
    #[clap(long, action)]
    report: bool,
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
        let out: Vec<Sequence> = seq_iter.collect();
        info!("{}",format!("{} sequences loaded from {}", out.len(), &seq_or_fname));
        return out;
    }

    return vec![Sequence::from_str(seq_name, seq_or_fname)];
}


pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();

    let mut reporters: Vec<Box<dyn AlignmentReporter>> = vec![];
    if args.pairwise { reporters.push(Box::new(PrintAsPairwise::new(80))); }
    if args.report { reporters.push(Box::new(SimilarityReport)); }
    if reporters.len() == 0 {  reporters.push(Box::new(SimilarityReport)); }

    let queries = get_sequences(&args.query, "query");

    if let Some(tmpl) = args.template {
        let templates = get_sequences(&tmpl, "template");
        align_all_pairs(&queries, &templates, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, false, &mut reporters);

    } else {
        align_all_pairs(&queries, &queries, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, true, &mut reporters);
    }
}