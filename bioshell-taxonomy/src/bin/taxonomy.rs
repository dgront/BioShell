use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use clap::Parser;
use bioshell_taxonomy::{Node, Rank, Taxonomy}; // Assuming your module is named `bioshell_taxonomy`
use std::path::PathBuf;
use log::info;
use env_logger;

use bioshell_io::{markdown_to_text};

const TAXONOMY_EXAMPLES: &str = include_str!("../documentation/taxonomy.md");
fn create_cookbook() -> String { format!("{}{}", "\x1B[4mCookbook:\x1B[0m\n", markdown_to_text(TAXONOMY_EXAMPLES)) }


/// Query NCBI Taxonomy from a local taxdump.tar.gz
#[derive(Parser, Debug)]
#[command(name = "taxonomy")]
#[clap(author, version, about, long_about = None, arg_required_else_help = true, after_long_help = create_cookbook())]
/// Provides taxonomy information
/// say taxonomy -h to see options
struct Cli {
    /// species name to lookup
    #[arg(short = 'n', long = "name")]
    name: Option<String>,

    /// file that lists species name to lookup
    #[arg(long)]
    names_file: Option<String>,

    /// taxonomic ID (taxid) to lookup
    #[arg(short = 't', long = "taxid")]
    taxid: Option<u32>,

    /// path to the taxonomy.dat or to taxdump.tar.gz file
    #[arg(short = 'p', long = "path", default_value = "./")]
    path: PathBuf,

    /// print also the selected nodes of the lineage for each species
    #[clap(short, long, short='l')]
    lineage: bool,

    /// print also the selected nodes of the lineage for each species
    #[clap(long)]
    lineage_full: bool,

    /// print the selected nodes of the lineage in columns
    #[clap(long)]
    lineage_table: bool,

    /// list known kingdoms
    #[clap(long)]
    kingdoms: bool,

    /// be more verbose and log program actions on the screen
    #[clap(short, long, short='v')]
    verbose: bool
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
        if args.verbose { env::set_var("RUST_LOG", "debug"); }
    }
    env_logger::init();

    let build_time = env!("BUILD_TIME");
    let git_commit_md5 = env!("GIT_COMMIT_MD5");

    info!("Build time: {}", build_time);
    info!("Git commit MD5 sum: {}", git_commit_md5);

    info!("Loading taxonomy from {}", args.path.display());
    let dump_file = format!("{}/taxdump.tar.gz", args.path.display());
    let taxonomy = Taxonomy::load_from_tar_gz(&dump_file)?;

    let mut nodes: Vec<&Node> = vec![];
    if let Some(name) = args.name {
        if let Some(taxid) = taxonomy.taxid(&name) {
            nodes.push(taxonomy.node(taxid).unwrap());
        } else {
            println!("Name '{}' not found in taxonomy.", name);
        }
    }

    if let Some(taxid) = args.taxid {
        if let Some(node) = taxonomy.node(taxid) {
            nodes.push(node);
        } else {
            println!("TaxID '{}' not found in taxonomy.", taxid);
        }
    }

    if let Some(fname) = args.names_file {
        let file = File::open(fname)?;
        let reader = BufReader::new(file);
        reader.lines().map(|line| line.unwrap()).for_each(|line| {
            if let Some(taxid) = taxonomy.taxid(&line) {
                if let Some(node) = taxonomy.node(taxid) {
                    nodes.push(node);
                }
            }
        });
    }
    if args.kingdoms {
        for n in taxonomy.nodes().filter(|ni|ni.rank==Rank::Kingdom) {
            let name = &n.name;
            let taxid = n.tax_id;
            print!("{} {} : ",n.tax_id, &name);
            for synonym in taxonomy.names(taxid) {
                if synonym != name { print!("{}; ", synonym); }
            }
            println!()
        }
    }
    for n in nodes {
        print!("{} {} : ", n.tax_id, n.name);
        for synonym in taxonomy.names(n.tax_id) {
            if synonym != &n.name { print!("{}; ", synonym); }
        }
        println!();

        if args.lineage {
            let lineage = taxonomy.lineage(n.tax_id);
            println!("Lineage:");
            for node in lineage.iter().filter(|n| n.rank != Rank::Unclassified && n.rank != Rank::Other) {
                println!("  {} ({:?})", node.name, node.rank);
            }
        }
        if args.lineage_table {
            let ranks = vec![Rank::Superkingdom,Rank::Kingdom, Rank::Phylum, Rank::Class, Rank::Order, Rank::Family, Rank::Genus, Rank::Species];
            for r in ranks {
                if let Some(rnk_node) = taxonomy.rank(n.tax_id, r) {
                    print!("{}\t", rnk_node.name);
                }
            }
            println!();
        }
        if args.lineage_full {
            let lineage = taxonomy.lineage(n.tax_id);
            println!("Lineage:");
            for node in lineage {
                println!("  {} ({:?})", node.name, node.rank);
            }
        }
    }

    Ok(())
}
