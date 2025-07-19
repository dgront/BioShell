use std::env;
use clap::{Parser};
#[allow(unused_imports)]
use log::{debug, info};

use bioshell_seq::sequence::filters::{AlwaysTrue, DescriptionContains, SequenceFilter};
use bioshell_seq::sequence::{clone_ungapped, count_identical, len_ungapped, ProfileColumnOrder,
            SequenceProfile, StockholmIterator};

use bioshell_io::{open_file, out_writer};
use bioshell_seq::msa::MSA;
use bioshell_seq::SequenceError;

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
    #[clap(short='i', long, action)]
    pairwise_identity: bool,
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

pub fn main() -> Result<(), SequenceError> {

    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();
    let args = Args::parse();

    let mut msa: MSA = MSA::default();
    // ---------- Read input MSA in the CLW (clustal-w) format
    if let Some(fname) = args.in_clw {
        let mut reader = open_file(&fname)?;
        let seq = StockholmIterator::from_stockholm_reader(&mut reader);
        msa = match MSA::from_sequences(seq) {
            Ok(msa) => msa,
            Err(error) => panic!("Incorrect sequence(s) found in MSA: {:?}", error),
        };
    }
    // ---------- Create the most probable sequence and print sequence identity between that sequence and any in the input set
    if args.max_sequence {
        let profile = SequenceProfile::new(ProfileColumnOrder::aa_standard_gapped(), &msa);
        let s = profile.most_probable_sequence();
        for si in msa.sequences() {
            let idnt = count_identical(&s, si).unwrap();
            let lenu =  len_ungapped(si);
            println!("{} {:4} / {:4} = {:.2}%", si.id(), idnt, lenu, idnt as f64 / lenu as f64 * 100.0);
        }
    }

    // ---------- Compute and print a sequence profile
    if args.profile {
        println!("{}",SequenceProfile::new(ProfileColumnOrder::aa_standard_gapped(), &msa));
    }

    if args.pairwise_identity {
        let max_seq_name_len: usize = msa.sequences().iter().map(|p| p.description().len()).max().unwrap_or(0).min(70);
        println!(
            "#{:^width$} {:^width$} {:^4} {:^4} {:^4} {:^6}",
            "i_seq_id", "j_seq_id", "idnt", "ilen", "jlen", "perc",
            width = max_seq_name_len
        );
        let mut min_id = 100.0_f64;
        let mut max_id = 0.0_f64;
        for si in msa.sequences() {
            let leni = len_ungapped(si);
            for sj in msa.sequences() {
                if si == sj { break; }
                let idnt = count_identical(&sj, si).unwrap();
                let lenj = len_ungapped(sj);
                let lenu = leni.min(lenj);
                let seq_id = idnt as f64 / lenu as f64 * 100.0;
                min_id = min_id.min(seq_id);
                max_id = max_id.max(seq_id);
                println!("{:width$} {:width$} {:4} {:4} {:4}  {:6.2}%", si.id(), sj.id(),
                         idnt, leni, lenj, seq_id, width = max_seq_name_len);
            }
        }
        println!("# min, max: {} {}", min_id, max_id);
    }

    // ---------- Convert the input MSA to FASTA format
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
    Ok(())
}