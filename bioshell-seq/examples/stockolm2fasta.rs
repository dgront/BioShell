use std::env;
use clap::{Parser};
use log::{info};

use bioshell_seq::sequence::{Sequence, remove_gaps_by_sequence, StockholmIterator};
use bioshell_io::open_file;

#[derive(Parser, Debug)]
#[clap(name = "stockholm2fasta")]
#[clap(about = "Converts a file in Stockholm format in FASTA", long_about = None)]
struct Args {
    /// input file in Stockholm format
    infile: String,
    /// remove gaps according to a sequence given by its ID
    #[clap(short='r')]
    if_remove: Option<String>
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();
    let fname= args.infile;

    let mut reader = open_file(&fname);

    let mut seq = StockholmIterator::from_stockholm_reader(&mut reader);

    if let Some(seq_id) = args.if_remove {
        info!("removing gapped column according to {}",seq_id);
        let mut ref_seq: Option<Sequence> = None;
        for seq in &seq {
            if seq.description()==seq_id { ref_seq = Option::from(seq.clone()); }
        }
        match ref_seq {
            None => { panic!("Can't find the reference sequence: {} in the input MSA",seq_id); }
            Some(the_ref_seq) => {
                remove_gaps_by_sequence(&the_ref_seq, &mut seq);
            }
        }
    }
    for s in &seq { println!("{}", s); }
    info!("{} sequences printed in FASTA format",seq.len());
}