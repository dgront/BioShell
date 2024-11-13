use std::env;
use clap::{Parser};
use log::{debug, info};
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
// use regex::Regex;

use bioshell_seq::sequence::{FastaIterator, count_residue_type, Sequence, ShorterThan, IsProtein, IsNucleic, SequenceFilter, LongerThan, ContainsX, LogicalNot, DescriptionContains};
use bioshell_io::{open_file, out_writer};
use bioshell_seq::SequenceError;

#[derive(Parser, Debug)]
#[clap(name = "filter_fasta")]
#[clap(about = "Filters sequences found in a FASTA file", long_about = None)]
struct Args {
    /// input file in FASTA format
    infile: String,
    /// remove sequences that are too short
    #[clap(short='l', long)]
    longer_than: Option<usize>,
    /// length of each line of the output sequence; use 0 to print the whole sequence on a single line
    #[clap(short='w', long, default_value_t = 80)]
    out_width: usize,
    /// remove sequences that are too long
    #[clap(short='s', long)]
    shorter_than: Option<usize>,
    /// remove sequences with more than given number of 'X' residues (unknowns)
    #[clap(short='x', long)]
    max_x: Option<usize>,
    /// print only valid protein sequences ('X' symbol is allowed)
    #[clap(short='p', long, action)]
    protein_only: bool,
    /// print only valid nucleotide sequences
    #[clap(short='n', long, action)]
    nucleotide_only: bool,
    /// print only unique sequences
    #[clap(short='u', long, action)]
    unique: bool,
    /// writes a file, where each row lists identical sequences; this works only when -u was used
    #[clap(long)]
    map_identical: Option<String>,
    /// batch retrieval by ID: print only these sequences whose ID is on the list provided as an input file
    #[clap(short='r', long)]
    retrieval_list: Option<String>,
    /// batch retrieval by subsequences: print a database sequence only if it contains aa subsequence from a given FASTA file; given strings must fully match a sequence's description to have it printed. For a substring match use the --description-contains option
    #[clap(short='m', long)]
    match_subsequences: Option<String>,
    /// print sequences whose description contains a certain substring, e.g. CYP1A or "DNA BINDING"
    #[clap(long)]
    description_contains: Option<String>,
    /// sort sequences by length
    #[clap(long, action)]
    sort: bool,

}

fn first_word_of_description(description: &str) -> String {
    let v: Vec<&str> = description.split_whitespace().collect();

    return String::from(v[0]);
}


pub fn main() -> Result<(), SequenceError> {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();
    let fname= args.infile;
    let out_width = args.out_width;

    let mut basic_filters: Vec<Box<dyn SequenceFilter>> = vec![];

    // ---------- Check is user wants to retrieve sequences by IDs
    let mut if_retrieve = false;
    let mut requested_ids: HashSet<String> = HashSet::new();
    if let Some(rfile) = args.retrieval_list {
        let file = File::open(rfile).expect("Can't read a list of seq-ids from a file");
        if_retrieve = true;
        let reader = BufReader::new(file);
        requested_ids = HashSet::from_iter(reader.lines().map(|l| l.unwrap()));
        info!("{} sequence IDs selected for retrieval", requested_ids.len());
    }

    // ---------- Check is user wants to retrieve sequences that are given in a query FASTA file
    let mut query_fasta: Vec<(String, String)> = vec![];
    let mut if_fasta_retrieve = false;
    if let Some(qfile) = args.match_subsequences {
        if_fasta_retrieve = true;
        let reader = open_file(&qfile)?;
        let seq_iter = FastaIterator::new(reader);
        for sequence in seq_iter {
            let id = String::from(sequence.id());
            query_fasta.push((sequence.to_string(out_width), id));
        }
    }


    let reader = open_file(&fname)?;
    let seq_iter = FastaIterator::new(reader);
    let mut cnt_all: usize = 0;
    let mut cnt_ok: usize = 0;
    // ---------- filter by length
    if let Some(min_length) = args.longer_than {
        basic_filters.push(Box::new(LongerThan { min_length }));
    }
    if let Some(max_length) = args.shorter_than {
        basic_filters.push(Box::new(ShorterThan { max_length }));
    }
    // ---------- filter by X
    if let Some(max_x) = args.max_x {
        basic_filters.push(Box::new(LogicalNot{ f: ContainsX{ min_x: max_x}}));
    }
    // ---------- filter by type: protein or nucleic acid
    if args.protein_only { basic_filters.push(Box::new(IsProtein)); }
    if args.nucleotide_only { basic_filters.push(Box::new(IsNucleic)); }
    // ---------- Check is user wants to filter sequences by a description substring
    if let Some(s) = args.description_contains {
        basic_filters.push(Box::new(DescriptionContains{ substring: s }))
    }

    let mut represented_sequences: HashMap<u64, Vec<String>> = HashMap::new();

    let mut output_seq: Option<Vec<Sequence>> = if args.sort { Some(vec![]) } else { None };

    let start = Instant::now();

    for sequence in seq_iter {
        cnt_all += 1;
        // ---------- check basic filters first
        let mut filters_ok = true;
        for filter in &basic_filters {
            if !filter.filter(&sequence) {
                filters_ok = false;
                break;
            }
        }
        if !filters_ok { continue }

        // ---------- keep only requested sequences
        if if_retrieve {
            let id = sequence.id();
            if !requested_ids.contains(id) { continue }
        }
        // ---------- keep only sequences found in the input query fasta
        if if_fasta_retrieve {
            let mut if_found = false;
            let sequence_as_str: String = sequence.to_string(out_width);
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
                    println!("{:width$}", s, width = out_width);
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
            let seq_id = first_word_of_description(&sequence.description());
            if represented_sequences.contains_key(&h) {
                represented_sequences.entry(h).or_default().push(seq_id);
                continue;
            }
            else {
                represented_sequences.insert(h, vec![seq_id]);
            }
        }
        cnt_ok += 1;
        if let Some(v) = &mut output_seq {
            v.push(sequence);
        } else {
            println!("{:width$}", sequence, width = out_width);
        }
    }

    if let Some(fname) = args.map_identical {
        if ! args.unique {
            panic!("The --map-identical option works only when --unique was also specified!");
        }
        let mut stream = out_writer(&fname, false);
        for (_hash, ids) in represented_sequences {
            let _ = stream.write(format!("{:}", ids[0]).as_bytes());
            for seq_id in &ids[1..] {
                let _ = stream.write(format!(" {:}", seq_id).as_bytes());
            }
            let _ = stream.write(b"\n");
        }
    }

    if let Some(v) = &mut output_seq {
        v.sort_by(|si, sj| {si.len().cmp(&sj.len())});
        for sequence in v { println!("{:width$}", sequence, width = out_width); }
    }

    info!("{} sequences processed in {:?}, {} of them printed in FASTA format",
        cnt_all, start.elapsed(), cnt_ok);

    Ok(())
}