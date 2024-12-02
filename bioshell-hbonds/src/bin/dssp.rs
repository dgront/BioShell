use std::env;
use clap::Parser;
use log::info;
use bioshell_hbonds::BackboneHBondMap;
use bioshell_pdb::Deposit;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None, arg_required_else_help = true)]
/// Command line tool to detect hydrogen bonds in a protein structure and to assign secondary structure
/// using the DSSP algorithm.
///
/// say dssp -h to see options
struct Args {
    /// input protein structure in either CIF or PDB format
    infile: String,
    /// print sequence and secondary structure in FASTA format
    #[clap(short, long, short='f')]
    out_fasta: bool,
    /// list all hydrogen bonds
    #[clap(long)]
    list: bool,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() {
    // ---------- BioShell app setup ----------
    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    if args.verbose {
        env::set_var("RUST_LOG", "debug");
    }
    env_logger::init();

    let build_time = env!("BUILD_TIME");
    let git_commit_md5 = env!("GIT_COMMIT_MD5");

    info!("Build time: {}", build_time);
    info!("Git commit MD5 sum: {}", git_commit_md5);

    // ---------- INPUT section ----------
    let deposit = Deposit::from_file(&args.infile).unwrap();
    let strctr= deposit.structure();

    // ---------- Detect H-bonds ----------
    let hbonds = BackboneHBondMap::new(&strctr);

    if args.list {
        for ((i,j), hbond) in hbonds.h_bonds() {
            println!("{:4} -> {:4}  {:.3}", i, j, hbond.dssp_energy());
        }
    }

}