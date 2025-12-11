#![doc(hidden)]

use std::env;
use clap::{Parser, ArgGroup};
use log::info;
use bioshell_io::{out_writer, markdown_to_text};
use bioshell_pdb::{Deposit, downlad_deposit_from_rcsb, EntityType, find_cif_file_name, find_pdb_file_name, PDBError, Structure};
use bioshell_pdb::pdb_atom_filters::{ByChain, ByEntity, InvertPredicate, IsBackbone, IsCA, IsHydrogen, IsNotWater, KeepNucleicAcid, KeepProtein, MatchAll, PdbAtomPredicate};
use bioshell_seq::chemical::ResidueTypeProperties;

mod deposit_info;

const PDB_TOOL_EXAMPLES: &str = include_str!("../documentation/pdb_tool.md");

fn create_cookbook() -> String { format!("{}{}", "\x1B[4mCookbook:\x1B[0m\n", markdown_to_text(PDB_TOOL_EXAMPLES)) }

#[derive(Parser, Debug)]
#[command(group(
    ArgGroup::new("input")
    .required(true)
    .args(&["infile", "pdb_id"])
))]
#[clap(author, version, about, long_about = None, arg_required_else_help = true, after_long_help = create_cookbook())]
/// Command line tool to operate on PDB files
/// say pdb_tool -h to see options
struct Args {
    /// input PDB file name
    #[clap(short, long, short='i')]
    infile: Option<String>,
    /// pdb_id of the input deposit
    #[clap(long)]
    pdb_id: Option<String>,
    /// path to a folder where a PDB / mmCIF file may be found
    #[clap(long, default_value="")]
    path: String,
    /// download the requested PDB deposit from the rcsb.org website
    #[clap(long)]
    online: bool,
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
    /// renumber residues of a structure
    #[clap(long, action)]
    renumber: bool,
    /// print a summary of the entities in the structure
    #[clap(long, action)]
    entities: bool,
    /// when printing entities, show also the sequence of an entity (only for polymer entities)
    #[clap(long, action)]
    entity_sequence: bool,
    /// keep only amino acid residues; all ligands and cofactors will be removed
    #[clap(long)]
    select_protein: bool,
    /// keep only nucleic acid residues; all proteins and ligands will be removed
    #[clap(long)]
    select_nucleic: bool,
    /// keep only selected chains
    #[clap(long)]
    select_chain: Option<String>,
    /// keep only alpha-carbon atoms
    #[clap(long, action)]
    select_ca: bool,
    /// keep only protein backbone atoms
    #[clap(long, action)]
    select_bb: bool,
    /// neglect water molecules while parsing the deposit
    #[clap(long, action)]
    skip_water: bool,
    /// neglect hydrogen atoms while parsing the deposit
    #[clap(long, action)]
    skip_hydrogens: bool,
    /// keep only selected entities
    #[clap(long, group = "select")]
    select_entity: Option<String>,
    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn filter<F: PdbAtomPredicate>(strctr: &Structure, filter: &F) -> Structure {

    let atoms_iter = strctr.atoms().iter().filter(|a| filter.check(a));
    return Structure::from_iterator(&strctr.id_code, atoms_iter.cloned());
}

fn print_info(deposit: &Deposit, tokens: &Vec<String>) {

    let tokens = deposit_info::get_deposit_info(deposit, tokens);
    for token in tokens {
        let val = token.1.replace([';', '\n'], " ");
        println!("{}: {}", token.0, val.trim());
    }
}

fn print_info_row(deposit: &Deposit, tokens: &Vec<String>) {

    let tokens = deposit_info::get_deposit_info(deposit, tokens);
    for token in tokens {
        let val = token.1.replace([';', '\n'], " ");
        print!("{}\t", val.trim());
    }
    println!();
}

fn write_pdb(strctr: &Structure, fname: &str) {
    let outstream = out_writer(fname, false);
    bioshell_pdb::write_pdb(strctr, outstream);
}

/// Print a list of all entities found in a structure.
///
/// This function lists only these entities that are present in the given structure; those present
/// in the deposit but not in the structure are not listed.
///
/// # Arguments
/// - substructure - structure to print entities for
/// - deposit - deposit containing the information about entities
/// - print_sequences - print also the sequence of an entity (only for polymer entities)
fn print_entities(substructure: &Structure, deposit: &Deposit, print_sequences: bool) {
    let entities = substructure.entity_ids();
    for (name, entity) in deposit.entities() {
        if !entities.contains(&name) { continue; }
        print!("{} {}", name, entity.description().trim());
        for chain in entity.chain_ids() {
            if substructure.chain_ids().contains(chain) { print!(" {}", chain); }
        }
        if let EntityType::Polymer(_poly_type) = entity.entity_type() {
            if let Some(db_ref) = entity.db_ref() {
                print!(" {}", &db_ref);
            }
            if print_sequences {
                let mut sequence : Vec<char> = vec![];
                for rt in  entity.entity_monomers() {
                    sequence.push(rt.parent_type.code1());
                }
                print!(" {}", Vec::from_iter(sequence).iter().collect::<String>());
            }
        }
        println!();
    }
}

fn load_deposit(args: &Args) -> Result<Deposit, PDBError> {

    // ---------- if a file name was given, load it
    if let Some(infile) = &args.infile {
        return Deposit::from_file(infile);
    }

    // ---------- if a PDB code was given, look for it online
    let pdb_id = args.pdb_id.as_ref().unwrap();
    if args.online {
        return downlad_deposit_from_rcsb(pdb_id);
    }

    // ---------- ... or look for a local file given the path
    let fname = find_cif_file_name(&pdb_id, &args.path)
        .or_else(|_| find_pdb_file_name(&pdb_id, &args.path))?;
    Deposit::from_file(fname)
}


fn main() -> Result<(), PDBError> {

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

    // ---------- INPUT section
    let deposit = load_deposit(&args)?;
    let mut strctr= deposit.structure()?;

    // ---------- FILTER section
    let mut multi_filter = MatchAll::new();
    if  args.select_protein {
        info!("Selecting only protein chains");
        multi_filter.add_predicate(Box::new(KeepProtein));
    }
    if args.select_nucleic {
        info!("Selecting only nucleic chains");
        multi_filter.add_predicate(Box::new(KeepNucleicAcid));
    }
    if let Some(chain_id) = args.select_chain {
        info!("Selecting only chain {}", &chain_id);
        multi_filter.add_predicate(Box::new(ByChain::new(&chain_id)));
    }
    if args.select_ca {
        info!("Selecting only alpha-carbon atoms");
        multi_filter.add_predicate(Box::new(IsCA));
    }
    if args.select_bb {
        info!("Selecting only proteins backbone atoms");
        multi_filter.add_predicate(Box::new(IsBackbone));
    }
    if args.skip_water {
        info!("Removing water molecules");
        multi_filter.add_predicate(Box::new(IsNotWater));
    }
    if args.skip_hydrogens {
        info!("Removing hydrogen atoms");
        multi_filter.add_predicate(Box::new(InvertPredicate::new(IsHydrogen)));
    }
    if let Some(entity_id) = args.select_entity {
        info!("Selecting entity: {}", entity_id);
        multi_filter.add_predicate(Box::new(ByEntity::new(&entity_id)));
    }
    if multi_filter.count_filters() > 0 {
            strctr = filter(&strctr, &multi_filter);
    }
    if args.renumber {
        strctr = strctr.renumbered_structure();
    }

    // ---------- OUTPUT section
    if args.out_fasta {
        for chain_id in &strctr.chain_ids() {
            println!("{}",strctr.sequence(chain_id));
        }
    }
    if args.entities {
        print_entities(&strctr, &deposit, args.entity_sequence);
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

    Ok(())
}

