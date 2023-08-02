use std::env;
use clap::{Parser};
use log::{debug, info};

use bioshell_seq::sequence::{clone_ungapped, count_identical, FastaIterator, len_ungapped, ProfileColumnOrder, remove_gaps, Sequence, SequenceProfile, StockholmIterator, AlwaysTrue, DescriptionContains, SequenceFilter};

use bioshell_io::{open_file, out_writer};
use bioshell_seq::msa::MSA;

#[derive(Parser, Debug)]
#[clap(name = "check_msa")]
#[clap(about = "Provides various statistics for a multiple sequence alignment", long_about = None)]
struct Args {
    /// input MSA file in the FASTA format
    #[clap(short='f', long)]
    in_fasta: Option<String>,
    /// input MSA file in the ClustalW or the Stockholm format
    #[clap(short='w', long)]
    in_clw: Option<String>,
    /// find the most probable sequence and compute its distance from each sequence of an MSA
    #[clap(short='a', long)]
    max_sequence: bool,
    /// compute sequence identity for all sequence pairs
    #[clap(short='i', long)]
    sequence_identity: bool,
    /// compute sequence profile
    #[clap(short='p', long)]
    profile: bool,
    /// remove gaps from all sequences
    #[clap(long)]
    remove_gaps: bool,
    /// select sequences which descriptions contain a given substring
    #[clap(short='s', long)]
    select: Option<String>,
    /// output MSA file in the FASTA format; use the --select option to save a part of the MSA
    #[clap(short='o', long)]
    out_fasta: Option<String>,
}

pub fn main() {

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args = Args::parse();

    let mut msa: MSA = MSA::default();
    if let Some(fname) = args.in_clw {
        let mut reader = open_file(&fname);
        let seq = StockholmIterator::from_stockholm_reader(&mut reader);
        msa = match MSA::from_sequences(seq) {
            Ok(msa) => msa,
            Err(error) => panic!("Incorrect sequence(s) found in MSA: {:?}", error),
        };
    }
    // ---------- Check is user wants to retrieve sequences by IDs
    if args.max_sequence {
        let profile = SequenceProfile::new(ProfileColumnOrder::aa_standard_gapped(), &msa);
        let s = profile.most_probable_sequence();
        for si in msa.sequences() {
            let idnt = count_identical(&s, si).unwrap();
            let lenu =  len_ungapped(si);
            println!("{} {:4} / {:4} = {:.2}%", si.id(), idnt, lenu, idnt as f64 / lenu as f64);
        }
    }

    // ---------- Compute and print a sequence profile
    if args.profile {
        println!("{}",SequenceProfile::new(ProfileColumnOrder::aa_standard_gapped(), &msa));
    }

    if let Some(fname) = args.out_fasta {
        let mut writer = out_writer(&fname, false);
        // ---------- Filter to print only selected sequences from the MSA
        let mut filter: Box<dyn SequenceFilter> = Box::new(AlwaysTrue);
        if let Some(substr) = args.select {
            filter = Box::new(DescriptionContains{ substring: substr });
        }

        for si in msa.sequences().iter().filter(|s| filter.filter(&s)) {
            if args.remove_gaps {
                writer.write(format!("{}", clone_ungapped(si)).as_bytes()).ok();
            } else {
                writer.write(format!("{}", si).as_bytes()).ok();
            }
        }
    }
}