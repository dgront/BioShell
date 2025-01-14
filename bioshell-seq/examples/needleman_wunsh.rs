use std::env;
use clap::{Parser};
#[allow(unused_imports)]
use log::{info};
use log::warn;
use bioshell_seq::sequence::{load_sequences};
use bioshell_seq::alignment::{PrintAsPairwise, SimilarityReport, align_all_pairs, MultiReporter, IdentityMatrixReporter};
use bioshell_seq::scoring::{SubstitutionMatrixList};
use bioshell_seq::SequenceError;

#[derive(Parser, Debug)]
#[clap(author, version,name = "needleman_wunsh")]
#[clap(allow_hyphen_values = true)]
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
    identity: bool,
    /// print sequence identity as a triangular matrix; this output is much more concise than the --identity output
    #[clap(long, action)]
    identity_matrix: bool,
    /// length of a sequence name to print; longer names will be trimmed that size
    #[clap(long, short='w', default_value = "20")]
    name_width: usize,
}


pub fn main() -> Result<(), SequenceError> {

    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }

    env_logger::init();
    let args = Args::parse();

    let name_width = args.name_width;
    let mut multireports = MultiReporter::new();
    if args.pairwise { multireports.add_reporter(Box::new(PrintAsPairwise::new(name_width, 80))); }
    if args.identity { multireports.add_reporter(Box::new(SimilarityReport::new(name_width))); }
    if args.identity_matrix {
        multireports.add_reporter(Box::new(IdentityMatrixReporter::new(name_width, "stdout")));
    }
    if multireports.count_reporters() == 0 {  multireports.add_reporter(Box::new(SimilarityReport::new(name_width))); }

    let queries = load_sequences(&args.query, "query")?;
    if queries.len() == 0 {
        warn!("No sequences found in the query set. Exiting.");
        return Ok(());
    }

    if let Some(tmpl) = args.template {
        let templates = load_sequences(&tmpl, "template")?;
        if templates.len() == 0 {
            warn!("No sequences found in the templates set. Exiting.");
            return Ok(());
        }
        align_all_pairs(&queries, &templates, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, false, &mut multireports);

    } else {
        align_all_pairs(&queries, &queries, SubstitutionMatrixList::BLOSUM62,
            args.open, args.extend, true, &mut multireports);
    }

    Ok(())
}