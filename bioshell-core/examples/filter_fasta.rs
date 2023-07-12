use std::env;
use clap::{Parser};
use log::{debug, info};
use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader};
use std::time::Instant;
use regex::Regex;

use bioshell_core::sequence::{FastaIterator, count_residue_type, Sequence};
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
    /// batch retrieval by ID: print only these sequences whose ID is on the list provided as an input file
    #[clap(short='r', long)]
    retrieval_list: Option<String>,
    /// batch retrieval by subsequences: print a database sequence only if it contains aa subsequence from a given FASTA file
    #[clap(short='m', long)]
    match_subsequences: Option<String>,
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();
    let fname= args.infile;

    // ---------- Check is user wants to retrieve sequences by IDs
    let mut if_retrieve = false;
    let mut requested_ids: HashSet<String> = HashSet::new();
    if let Some(rfile) = args.retrieval_list {
        let file = File::open(rfile).expect("Can't read a list of seq-ids from a file");
        if_retrieve = true;
        let reader = BufReader::new(file);
        requested_ids = HashSet::from_iter(reader.lines().map(|l| l.unwrap()));
        info!("{} sequences selected for retrieval", requested_ids.len());
    }

    // ---------- Check is user wants to retrieve sequences that are given in a query FASTA file
    let mut query_fasta: Vec<(String, String)> = vec![];
    let mut if_fasta_retrieve = false;
    if let Some(qfile) = args.match_subsequences {
        if_fasta_retrieve = true;
        let reader = open_file(&qfile);
        let seq_iter = FastaIterator::new(reader);
        for sequence in seq_iter {
            let id = String::from(sequence.id());
            query_fasta.push((sequence.to_string(), id));
        }
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
        // ---------- filter by length
        if sequence.len() < min_len || sequence.len() > max_len { continue }
        // ---------- remove too many X's
        if max_x > 0 && count_residue_type(&sequence, 'X') > max_x { continue }
        // ---------- keep only requested sequences
        if if_retrieve {
            let id = sequence.id();
            if !requested_ids.contains(id) { continue }
        }
        // ---------- keep only sequences found in the input query fasta
        if if_fasta_retrieve {
            let mut if_found = false;
            let sequence_as_str: String = sequence.to_string();
            for (seq, id) in &query_fasta {
                let it: Vec<_>  = sequence_as_str.match_indices(seq).collect();
                if it.len() > 0 {
                    if_found = true;
                    debug!("Found {} in {}",id, sequence.id());
                }
                for (from, hit) in it {
                    let perc = (hit.len() as f64) / (sequence.len() as f64) * 100.0;
                    let new_id = format!("{} ({} {}:{}, {:6.2}%)",
                                         sequence.description(), id, from, from + hit.len(), perc);
                    let s = Sequence::new(&new_id, &sequence_as_str);
                    println!("{}", s);
                }
            }
            if !if_found { debug!("Nothing found for {}",sequence.id()); }
            if cnt_all % 100 == 0 { debug!("Processed {} sequences",cnt_all); }
            continue;   // --- don't go below, since the sequence is already printed; instead start with the next sequence
        }
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