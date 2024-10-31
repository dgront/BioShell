use std::env;
use clap::{Parser};
use log::info;
use bioshell_io::{out_writer};
use bioshell_pdb::{Deposit, Structure};
use bioshell_pdb::pdb_atom_filters::{ByChain, ByEntity, IsCA, MatchAll, PdbAtomPredicate};

mod deposit_info;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None, arg_required_else_help = true)]
/// Command line tool to operate on PDB files
/// say pdb_tool -h to see options
struct Args {
    /// input PDB file name
    #[clap(short, long, short='i', required=true)]
    infile: String,
    /// print FASTA sequence for every chain in each input file
    #[clap(short, long, short='f')]
    out_fasta: bool,
    /// print selected structure in PDB format
    #[clap(short, long)]
    out_pdb: Option<String>,
    /// print secondary structure for every chain in each input file
    #[clap(long)]
    out_secondary: bool,
    /// print basic information about a given structure
    ///
    /// the following parameters can be used to select which information to print:
    /// - 'id' - deposit ID code
    /// - 'keywords' - list of keywords
    /// - 'title' - deposit title
    /// - 'resolution' - resolution
    /// - 'methods' - list of methods
    /// - 'ligands' - list of ligands
    /// - 'rfactor' - deposit keywords
    /// - 'spacegroup' - resolution
    /// - 'unitcell' - list of methods
    /// - 'classification' - list of ligands
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 0..)]
    info: Option<Vec<String>>,
    /// print basic information about a given structure in the tsv format.
    ///
    /// The option accepts the same keys as the 'info' option
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 0..)]
    info_table: Option<Vec<String>>,
    /// break FASTA lines when longer that given cutoff
    #[clap(long, default_value="80")]
    out_fasta_width: usize,
    /// keep only selected chains
    #[clap(long)]
    select_chain: Option<String>,
    /// keep only alpha-carbon atoms
    #[clap(long)]
    select_ca: bool,
    /// keep only selected entities
    #[clap(long, group = "select")]
    select_entity: Option<String>,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn filter<F: PdbAtomPredicate>(strctr: &Structure, filter: &F) -> Structure {

    let atoms_iter = strctr.atoms().iter().filter(|a| filter.check(a));
    return Structure::from_iterator(&strctr.id_code, atoms_iter);
}

fn print_info(deposit: &Deposit, tokens: &Vec<String>) {

    let tokens = deposit_info::get_deposit_info(deposit, tokens);
    for token in tokens {
        println!("{}: {}", token.0, token.1);
    }
}

fn print_info_row(deposit: &Deposit, tokens: &Vec<String>) {

    let tokens = deposit_info::get_deposit_info(deposit, tokens);
    for token in tokens {
        print!("{}\t", token.1);
    }
    println!();
}

fn write_pdb(strctr: &Structure, fname: &str) {
    let mut outstream = out_writer(fname, false);
    for a in strctr.atoms() { write!(outstream, "{}\n", a).unwrap(); }
    outstream.flush().unwrap();
}

fn main() {
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

    // ---------- INPUT section
    let deposit = Deposit::from_file(&args.infile).unwrap();
    let mut strctr= deposit.structure();

    // ---------- FILTER section
    let mut multi_filter = MatchAll::new();
    if let Some(chain_id) = args.select_chain {
        info!("Selecting only chain {}", &chain_id);
        multi_filter.add_predicate(Box::new(ByChain::new(&chain_id)));
    }
    if args.select_ca {
        info!("Selecting only alpha-carbon atoms");
        multi_filter.add_predicate(Box::new(IsCA));
    }
    if let Some(entity_id) = args.select_entity {
        info!("Selecting entity: {}", entity_id);
        multi_filter.add_predicate(Box::new(ByEntity::new(&entity_id)));
    }
    if multi_filter.count_filters() > 0 {
        strctr = filter(&strctr, &multi_filter);
    }

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
    if let Some(tokens) = args.info {
        print_info(&deposit, &tokens);
    }
    if let Some(tokens) = args.info_table {
        print_info_row(&deposit, &tokens);
    }

    if let Some(out_fname) = args.out_pdb {
        write_pdb(&strctr, &out_fname);
    }
}

