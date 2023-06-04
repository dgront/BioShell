use std::env;
use clap::{Parser};
use log::{info};
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader};
use std::time::Instant;

use bioshell_core::sequence::{FastaIterator, count_residue_type};
use bioshell_core::utils::open_file;

#[derive(Parser, Debug)]
#[clap(name = "filter_fasta")]
#[clap(about = "Filters sequences found in a FASTA file", long_about = None)]
struct Args {
    /// input file in FASTA format
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
    /// print only unique sequences
    #[clap(short='u', long)]
    unique: bool,
    /// batch retrieval: print only these sequences whose ID is on the list provided as an input file
    #[clap(short='r', long)]
    retrieval_list: Option<String>,
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();
    let fname= args.infile;
    let mut if_retrieve = false;
    let mut requested_ids: HashSet<String> = HashSet::new();
    if let Some(rfile) = args.retrieval_list {
        let file = File::open(rfile).expect("Can't read a list of seq-ids from a file");
        if_retrieve = true;
        let reader = BufReader::new(file);
        requested_ids = HashSet::from_iter(reader.lines().map(|l| l.unwrap()));
        info!("{} sequences selected for retrieval", requested_ids.len());
    }
    let reader = open_file(&fname);
    let seq_iter = FastaIterator::new(reader);
    let mut cnt_all: usize = 0;
    let mut cnt_ok: usize = 0;
    let min_len = match args.longer_than { None => 0, Some(l) => l };
    let max_len = match args.shorter_than { None => 1000000, Some(l) => l };
    let max_x = match args.max_x { None => 0, Some(l) => l };

    let mut observed_sequences:HashSet<u64> = HashSet::new();

    let start = Instant::now();
    for sequence in seq_iter {
        cnt_all += 1;
        // ---------- keep only requested sequences
        if if_retrieve {
            let id = sequence.id();
            if !requested_ids.contains(id) { continue }
        }
        // ---------- filter by length
        if sequence.len() < min_len || sequence.len() > max_len { continue }
        // ---------- remove too many X's
        if max_x > 0 && count_residue_type(&sequence, 'X') > max_x { continue }
        // ---------- remove redundant sequences
        if args.unique {
            let mut hasher = DefaultHasher::new();
            sequence.hash(&mut hasher);
            let h = hasher.finish();
            if observed_sequences.contains(&h) { continue}
            else {observed_sequences.insert(h);}
        }
        cnt_ok += 1;
        println!("{}", sequence);
    }

    info!("{} sequences processed in {:?}, {} of them printed in FASTA format",
        cnt_all, start.elapsed(), cnt_ok);
}