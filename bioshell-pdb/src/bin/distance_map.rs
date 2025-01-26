use std::env;
use clap::Parser;
use bioshell_pdb::{Deposit, PdbAtom, Structure};
use bioshell_pdb::calc::distance;
use bioshell_pdb::pdb_atom_filters::{AlwaysPass, ByChain, IsBackbone, IsCA, IsCB, KeepProtein, MatchAll, PdbAtomPredicate};
use log::info;

fn describe_atom(a: &PdbAtom) -> String {
    format!("{} {} {:4} {}", a.chain_id, a.res_name, a.res_seq, a.name)
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None, arg_required_else_help = true)]
/// Simple program to compute a matrix of inter-atomic distances for a given PDB file
/// say distance_map -h to see options
struct Args {
    /// input PDB file name
    #[clap(short, long, short='i', required=true)]
    infile: String,
    /// use only C-alpha atoms for calculations
    #[clap(long)]
    ca: bool,
    /// use only backbone atoms for calculations
    #[clap(long)]
    bb: bool,
    /// use only C-beta atoms for calculations
    #[clap(long)]
    cb: bool,
    /// keep only amino acid residues; all ligands and cofactors will be removed
    #[clap(long)]
    select_protein: bool,
    /// keep only selected chain
    #[clap(long)]
    select_chain: Option<String>,
    /// print contact map only
    #[clap(long)]
    cmap: Option<f64>,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn print_distance_map(structure: &Structure) {

    for ai in structure.atoms().iter() {
        for aj in structure.atoms().iter() {
            if ai.res_seq == aj.res_seq && ai.i_code==aj.i_code && ai.chain_id==aj.chain_id { break }
            println!("{} {} : {:.3}", describe_atom(ai), describe_atom(aj), distance(ai, aj));
        }
    }
}

fn print_contacts(structure: &Structure, cutoff: f64) {

    for ai in structure.atoms().iter() {
        for aj in structure.atoms().iter() {
            if ai.res_seq == aj.res_seq && ai.i_code==aj.i_code && ai.chain_id==aj.chain_id { break }
            if distance(ai, aj) <= cutoff {println!("{} {}", describe_atom(ai), describe_atom(aj)); }
        }
    }
}

fn main() {
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
    // ---------- FILTER section
    let mut multi_filter = MatchAll::new();
    let mut atom_selector: Box<dyn PdbAtomPredicate> = Box::new(AlwaysPass);
    if args.bb { atom_selector = Box::new(IsBackbone); }
    if args.ca { atom_selector = Box::new(IsCA); }
    if args.cb {atom_selector = Box::new(IsCB); }
    multi_filter.add_predicate(atom_selector);
    if  args.select_protein {
        info!("Selecting only protein chains");
        multi_filter.add_predicate(Box::new(KeepProtein));
    }
    if let Some(chain_id) = args.select_chain {
        info!("Selecting only chain {}", &chain_id);
        multi_filter.add_predicate(Box::new(ByChain::new(&chain_id)));
    }

    let deposit = Deposit::from_file(&args.infile).unwrap();
    let strctr = deposit.structure().unwrap();

    let atoms_iter = strctr.atoms().iter().filter(|a| multi_filter.check(a));
    let strctr = Structure::from_iterator(&strctr.id_code, atoms_iter);

    if let Some(cutoff) = args.cmap {
        print_contacts(&strctr, cutoff);
    } else {
        print_distance_map(&strctr);
    }
}
