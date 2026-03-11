use std::env;
use clap::{Parser};
#[allow(unused_imports)]
use log::{debug, info};

use bioshell_seq::sequence::filters::{DescriptionContains, SequenceFilter};
use bioshell_seq::sequence::{count_identical, FastaIterator, len_ungapped, ProfileColumnOrder, remove_gaps_by_sequence, Sequence, SequenceProfile, trim_by_sequence};

use bioshell_io::{open_file, out_writer};
use bioshell_seq::msa::{MSA, StockholmMSA};
use bioshell_seq::SequenceError;

#[derive(Parser, Debug)]
#[clap(name = "check_msa")]
#[clap(about = "Provides various statistics for a multiple sequence alignment", long_about = None)]
struct Args {
    /// input MSA file in the FASTA format
    #[clap(short='f', long)]
    in_fasta: Option<String>,
    /// input MSA file in the Stockholm format
    #[clap(short='s', long)]
    in_sto: Option<String>,
    /// prints info about number of sequences in a given MSA and their lengths
    #[clap(long, action)]
    info: bool,
    /// find the most probable sequence and compute its distance from each sequence of an MSA
    #[clap(short='a', long)]
    max_sequence: bool,
    /// compute sequence identity for all sequence pairs
    #[clap(short='i', long, action)]
    pairwise_identity: bool,
    /// compute sequence profile
    #[clap(short='p', long)]
    profile: bool,
    /// select sequences which descriptions contain a given substring
    #[clap(short='d', long)]
    desc_has: Option<String>,
    /// remove all the letters from all sequences that are aligned to a gap in the reference sequence
    ///
    /// This option creates an alignment where the reference sequence (pointed out by its ID) will
    /// have no gaps and all other sequences may be shortened
    #[clap(long)]
    trim_by_gaps: Option<String>,
    /// remove all the letters from both ends of all sequences that are aligned to an N-terminal or C-terminal gap in the reference sequence
    #[clap(long)]
    trim_by_ends: Option<String>,
    /// output MSA file in the FASTA format
    #[clap(short='o', long)]
    out_fasta: Option<String>,
    /// output MSA file in the Stockholm format
    #[clap(long)]
    out_sto: Option<String>,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool,
}

fn print_info(msa: &MSA) {

    let n_seq = msa.n_seq();
    println!("n_sequences: {:4}", n_seq);
    if n_seq == 0 { return; }
    let mut lengths: Vec<usize> = vec![];
    lengths.reserve(n_seq);
    for seq in msa.sequences() {
        lengths.push(len_ungapped(seq));
    }
    lengths.sort();
    println!("n_columns:   {:4}", msa.sequences()[0].len());
    println!("shortest:    {:4}", lengths[0]);
    println!("longest:     {:4}", lengths[n_seq-1]);
    println!("median:      {:4}", lengths[n_seq/2]);
}


pub fn main() -> Result<(), SequenceError> {

    let args = Args::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    let build_time = env!("BUILD_TIME");
    let git_commit_md5 = env!("GIT_COMMIT_MD5");

    info!("Build time: {}", build_time);
    info!("Git commit MD5 sum: {}", git_commit_md5);


    // ---------- Create an empty MSA
    let mut msa: StockholmMSA = MSA::default().into();

    // ---------- Read input MSA in the Stockholm (.sto) format
    if let Some(fname) = args.in_sto {
        let mut reader = open_file(&fname)?;
        msa = StockholmMSA::from_stockholm_reader(&mut reader)?;
    }
    // ---------- Read input MSA in the .fasta format
    if let Some(fname) = args.in_fasta {
        let reader = open_file(&fname)?;
        let seq: Vec<Sequence> = FastaIterator::new(reader).collect();
        msa = match MSA::from_sequences(seq) {
            Ok(msa) => msa.into(),
            Err(error) => panic!("Incorrect sequence(s) found in MSA: {:?}", error),
        };
    }

    // ---------- Print basic information about the given MSA
    if args.info {
        print_info(&msa);
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

    // ---------- Trim the alignment by the reference sequence
    if let Some(seq_id) = args.trim_by_gaps {
        let ref_seq = msa.by_id(&seq_id).ok_or(SequenceError::InvalidSequenceID{ seq_id: seq_id.to_string() })?;
        info!("removing gapped columns according to {}",seq_id);
        let mut filtered_seq = remove_gaps_by_sequence(ref_seq, msa.sequences());
        msa.clear();
        for seq in filtered_seq.into_iter() {
            msa.add_sequence(seq)?;
        }
    }

    // ---------- Trim the alignment by the reference sequence
    if let Some(seq_id) = args.trim_by_ends {
        let ref_seq = msa.by_id(&seq_id).ok_or(SequenceError::InvalidSequenceID{ seq_id: seq_id.to_string() })?;
        info!("trimming the MSA according to {}",seq_id);
        let mut filtered_seq = trim_by_sequence(ref_seq, msa.sequences())?;
        msa.clear();
        for seq in filtered_seq.into_iter() {
            msa.add_sequence(seq)?;
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

    if let Some(substr) = args.desc_has {
        let filter = Box::new(DescriptionContains{ substring: substr });
        let seq: Vec<Sequence> = msa.sequences().iter().filter(|s| filter.filter(s)).cloned().collect();
        msa = MSA::from_sequences(seq)?.into();
    }

    // ---------- Convert the input MSA to Stockholm format
    if let Some(fname) = args.out_sto {
        let mut writer = out_writer(&fname, false);
        write!(writer, "{}", msa)?;
    }

    // ---------- Convert the input MSA to FASTA format
    if let Some(fname) = args.out_fasta {
        let mut writer = out_writer(&fname, false);
        write!(writer, "{}", msa.msa())?;
    }
    Ok(())
}