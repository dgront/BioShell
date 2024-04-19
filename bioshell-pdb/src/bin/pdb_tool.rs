use std::env;
use clap::{Parser};
use bioshell_pdb::{load_pdb_file, Structure};
use bioshell_pdb::pdb_atom_filters::{ByChain, IsCA, PdbAtomPredicate};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Command line tool to operate on PDB files
/// say pdb_tool -h to see options
struct Args {
    /// input PDB file name
    #[clap(short, long, short='i', required=true)]
    infile: String,
    /// print FASTA sequence for every chain in each input file
    #[clap(short, long, short='f')]
    out_fasta: bool,
    /// print secondary structure for every chain in each input file
    #[clap(long)]
    out_secondary: bool,
    /// print basic information about a given structure
    #[clap(long)]
    info: bool,
    /// break FASTA lines when longer that given cutoff
    #[clap(long, default_value="80")]
    out_fasta_width: usize,
    /// keep only selected chains
    #[clap(long, group = "select")]
    select_chain: Option<String>,
    /// keep only alpha-carbond
    #[clap(long, group = "select")]
    select_ca: bool,
}

fn filter(strctr: &mut Structure, filters: &Vec<Box<dyn PdbAtomPredicate>>) {
    for f in filters {
        strctr.atoms_mut().retain(|a| !f.check(&a));
    }
}

fn main() {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // ---------- INPUT section
    let mut strctr = load_pdb_file(&args.infile).unwrap();

    // ---------- FILTER section
    let mut filters: Vec<Box<dyn PdbAtomPredicate>> = vec![];
    if let Some(chain_id) = args.select_chain {
        filters.push(Box::new(ByChain::new(&chain_id)));
    }
    if args.select_ca { filters.push(Box::new(IsCA)); }
    filter(&mut strctr, &filters);

    // ---------- OUTPUT section
    if args.out_fasta {
        for chain_id in &strctr.chain_ids() {
            println!("{}",strctr.sequence(chain_id));
        }
    }
    if args.out_secondary {
        for chain_id in &strctr.chain_ids() {
            let seq = strctr.sequence(chain_id);
            println!("> {}\n{}", seq.description(), strctr.secondary(chain_id).to_string());
        }
    }
    if args.info {
        println!("methods: {:?}",strctr.methods);
    }
}

