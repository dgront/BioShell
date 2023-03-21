use clap::Parser;
use log::info;
use std::env;

use bioshell_core::sequence::{
    from_file, from_stockholm_reader, remove_gaps_by_sequence, Sequence,
};

#[derive(Parser, Debug)]
#[clap(name = "stockholm2fasta")]
#[clap(about = "Converts a file in Stockholm format in FASTA", long_about = None)]
struct Args {
    /// input file in Stockholm format
    infile: String,
    #[clap(short = 'r')]
    if_remove: Option<String>,
}

pub fn main() {
    if env::var("RUST_LOG").is_err() {
        env::set_var("RUST_LOG", "info")
    }
    let args = Args::parse();
    let mut seq: Vec<Sequence> = from_file(&args.infile, from_stockholm_reader);
    if let Some(seq_id) = args.if_remove {
        info!("removing gapped column according to {}", seq_id);
        let mut ref_seq: Option<Sequence> = None;
        for seq in &seq {
            if seq.id() == seq_id {
                ref_seq = Option::from(seq.clone());
            }
        }
        match ref_seq {
            None => {
                panic!(
                    "Can't find the reference sequence: {} in the input MSA",
                    seq_id
                );
            }
            Some(the_ref_seq) => {
                remove_gaps_by_sequence(&the_ref_seq, &mut seq);
            }
        }
    }
    for s in &seq {
        println!("{}", s);
    }
    info!("{} sequences printed in FASTA format", seq.len());
}
