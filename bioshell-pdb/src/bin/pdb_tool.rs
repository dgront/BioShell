use std::env;
use clap::{Parser};
use log::info;
use bioshell_cif::is_cif_file;
use bioshell_io::{open_file, out_writer};
use bioshell_pdb::{is_pdb_file, load_cif_file, load_pdb_file, Structure};
use bioshell_pdb::pdb_atom_filters::{ByChain, ByEntity, IsCA, MatchAll, PdbAtomPredicate};

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
    #[clap(long)]
    info: bool,
    /// print basic information about a given structure in the tsv format
    #[clap(long)]
    info_table: bool,
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

fn print_info(strctr: &Structure) {
    println!("id_code: {:?}",strctr.id_code);
    println!("methods: {:?}",strctr.methods);
    if let Some(class) = &strctr.classification { println!("classification: {}", class); }
    if let Some(title) = &strctr.title { println!("title: {}", title); }
    if let Some(res) = strctr.resolution { println!("resolution: {}", res); }
    if let Some(r_fact) = strctr.r_factor { println!("r_factor: {}", r_fact); }
    if let Some(r_free) = strctr.r_free { println!("r_free: {}", r_free); }
    if let Some(unit_cell) = &strctr.unit_cell { println!("space group: {}", unit_cell.space_group); }
    if strctr.keywords.len() > 0 { println!("keywords: {}", strctr.keywords.join(", ")); }
    println!("models: {}", strctr.count_models());
}

fn print_info_row(strctr: &Structure) {
    print!("{:?}\t",strctr.id_code);
    print!("{:?}\t",strctr.methods);
    if let Some(class) = &strctr.classification { 
        print!("{}\t", class); 
    } else { print!("\t"); }
    if let Some(title) = &strctr.title { 
        print!("{}\t", title.replace("\n", " ")); 
    } else { print!("\t"); }
    if let Some(res) = strctr.resolution { 
        print!("{}\t", res); 
    } else { print!("\t"); }
    if let Some(r_fact) = strctr.r_factor { 
        print!("{}\t", r_fact); 
    } else { print!("\t"); }
    if let Some(r_free) = strctr.r_free { 
        print!("{}\t", r_free); 
    } else { print!("\t"); }
    if let Some(unit_cell) = &strctr.unit_cell { 
        print!("{}\t", unit_cell.space_group); 
    } else { print!("\t"); }
    print!("{}\n", strctr.count_models());
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
    let mut strctr: Structure;
    if is_pdb_file(&args.infile).is_ok_and(|f| f) {
        strctr = load_pdb_file(&args.infile).unwrap();
    } else if is_cif_file(&args.infile).is_ok_and(|f| f) {
        strctr = load_cif_file(&args.infile).unwrap();
    } else {
        panic!("Can't recognize the format of the input file: {}", &args.infile);
    }

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
    if args.info { print_info(&strctr); }
    if args.info_table { print_info_row(&strctr); }

    if let Some(out_fname) = args.out_pdb {
        write_pdb(&strctr, &out_fname);
    }
}

