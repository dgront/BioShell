use std::env;
use clap::{Parser};
use log::{info};

use bioshell_core::sequence::{Sequence, FastaIterator, count_residue_type};
use bioshell_core::utils::open_file;

#[derive(Parser, Debug)]
#[clap(name = "stockholm2fasta")]
#[clap(about = "Converts a file in Stockholm format in FASTA", long_about = None)]
struct Args {
    /// input file in Stockholm format
    infile: String,
    /// remove sequences that are too short
    #[clap(short='l', long)]
    longer_than: Option<usize>,
    /// remove sequences that are too long
    #[clap(short='s', long)]
    shorter_than: Option<usize>,
    /// remove sequences with more than given number of 'X' residues (unknowns)
    #[clap(short='x', long)]
    max_x: Option<usize>,
    /// print unique sequences
    #[clap(short='u', long)]
    unique: bool,
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    let args = Args::parse();
    let fname= args.infile;

    let reader = open_file(&fname);
    let seq_iter = FastaIterator::new(reader);
    let mut cnt: usize = 0;
    let min_len = match args.longer_than { None => 1000000000, Some(l) => l };
    let max_len = match args.shorter_than { None => 0, Some(l) => l };
    let max_x = match args.max_x { None => 0, Some(l) => l };

    for sequence in seq_iter {
        // ---------- filter by length
        if sequence.len() < min_len || sequence.len() > max_len { continue }
        // ---------- remove too many X's
        if max_x > 0 && count_residue_type(&sequence, 'X') > max_x { continue }
        cnt += 1;
        println!("{}", sequence);
    }
    info!("{} sequences printed in FASTA format",cnt);
}