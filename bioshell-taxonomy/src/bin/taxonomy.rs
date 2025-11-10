use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use clap::Parser;
use bioshell_taxonomy::{Node, Rank, Taxonomy, TaxonomyMatcher}; // Assuming your module is named `bioshell_taxonomy`
use std::path::PathBuf;
use log::info;
use env_logger;

use bioshell_io::{markdown_to_text, open_file};
use bioshell_seq::sequence::{FastaIterator, Sequence, SequenceReporter, WriteFasta};

const TAXONOMY_EXAMPLES: &str = include_str!("../documentation/taxonomy_app.md");
fn create_cookbook() -> String { format!("{}{}", "\x1B[4mCookbook:\x1B[0m\n", markdown_to_text(TAXONOMY_EXAMPLES)) }


#[derive(Parser, Debug)]
#[command(name = "taxonomy")]
#[clap(author, version, about, long_about = None, arg_required_else_help = true, after_long_help = create_cookbook())]
/// Provides taxonomy information based on  NCBI taxonomy data loaded from a local taxdump.tar.gz file
/// say taxonomy -h to see options
struct Cli {
    /// species name to lookup
    #[arg(short = 'n', long = "name")]
    name: Option<String>,

    /// detect species in a sequence description
    #[arg(short = 'd', long)]
    detect: Option<String>,

    /// detect species in every sequence description in a given .fasta file
    #[arg(short = 'f', long)]
    detect_fasta: Option<String>,

    /// file that lists species name to lookup
    #[arg(long)]
    names_file: Option<String>,

    /// taxonomic ID (taxid) to lookup
    #[arg(short = 't', long = "taxid")]
    taxid: Option<u32>,

    /// file that lists taxonomic ID (taxid) to lookup
    #[arg(long)]
    taxid_file: Option<String>,

    /// path to the taxonomy.dat or to taxdump.tar.gz file
    #[arg(short = 'p', long = "path", default_value = "./")]
    path: PathBuf,

    /// download the most recent taxdump.tar.gz file from the NCBI website
    #[clap(long)]
    download: bool,

    /// string used as a separator for the text fields in the output; the default is ' : ' (space colon space)
    #[arg(long = "separator")]
    separator: Option<String>,

    /// print also the selected nodes of the lineage for each species
    #[clap(short, long, short='l')]
    lineage: bool,

    /// print also the selected nodes of the lineage for each species
    #[clap(long)]
    lineage_full: bool,

    /// print the selected nodes of the lineage in columns
    #[clap(long)]
    lineage_table: bool,

    /// print the kingdom each taxid / name belongs to
    #[clap(long)]
    kingdom: bool,

    /// print the domain (superkingdom) each taxid / name  belongs to
    #[clap(long)]
    domain: bool,

    /// print the rank for each requested each taxid / name
    #[clap(long)]
    rank: bool,

    /// print the results as JSON rather than as a free text
    #[clap(long)]
    json: bool,

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

    if args.download {
        Taxonomy::download_from_ncbi("taxdump.tar.gz")?;
    }

    info!("Loading taxonomy from {}", args.path.display());
    let dump_file = format!("{}/taxdump.tar.gz", args.path.display());
    let taxonomy = Taxonomy::load_from_tar_gz(&dump_file)?;

    let separator: &str = if let Some(ref sep) = args.separator { sep } else { " : " };
        
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
    if let Some(fname) = args.taxid_file {
        let file = File::open(fname)?;
        let reader = BufReader::new(file);
        for taxid in reader.lines().map(|line| line.unwrap().parse::<u32>().unwrap()) {
            if let Some(node) = taxonomy.node(taxid) {
                nodes.push(node);
            }
        }
    }

    for n in &nodes {
        if args.json {
            if !args.lineage {
                println!("{}", serde_json::to_string_pretty(&n).expect("Can't serialize a Node struct!"));
            }
        } else {
            print!("{}{}{}{} ", n.tax_id, separator, n.name, separator);
        }
        if args.rank { print!("{:?}{}", n.rank, separator); }
        if args.kingdom {
            if let Some(kingdom) = taxonomy.rank(n.tax_id, Rank::Kingdom) { print!("{}{}", kingdom.name, separator); }
        }
        if args.domain {
            if let Some(domain) = taxonomy.rank(n.tax_id, Rank::Superkingdom) { print!("{}{}", domain.name, separator); }
        }
        if !args.json {
//            print!("{}", separator);

            for synonym in taxonomy.names(n.tax_id) {
                if synonym != &n.name { print!("{}{} ", synonym, separator); }
            }
            println!();
        }
        if args.lineage {
            let lineage = taxonomy.lineage(n.tax_id);
            if args.json {
                println!("{}",serde_json::to_string_pretty(&lineage).expect("Can't serialize a Node struct!"));
            } else {
                println!("Lineage:");
                for node in lineage.iter().filter(|n| n.rank != Rank::Unclassified && n.rank != Rank::Other) {
                    println!("  {} ({:?})", node.name, node.rank);
                }
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

    if args.detect.is_none() && args.detect_fasta.is_none() {
        return Ok(());
    }

    let matcher: TaxonomyMatcher = TaxonomyMatcher::from_taxonomy(taxonomy)?;
    let taxonomy = matcher.taxonomy();

    if let Some(desc) = args.detect {
        let taxid = matcher.find(&desc);
        if let Some(taxid) = taxid {
            let node = taxonomy.node(taxid).unwrap();
            println!("TaxId={} {}", node.tax_id, node.name);
        }
    }

    if let Some(fasta) = args.detect_fasta {
        let mut reporter = WriteFasta::new(None, 0, false);
        let sequences = FastaIterator::new(open_file(fasta)?);
        for seq in sequences {
            let desc = seq.description();
            if let Some(taxid) = matcher.find(&desc) {
                let node = taxonomy.node(taxid).unwrap();
                let new_desc = format!("{} TaxId={}[{}]", desc, node.tax_id, node.name);
                let s = Sequence::new(&new_desc, std::str::from_utf8(seq.as_u8())?);
                reporter.report(&s)?;
            } else {
                reporter.report(&seq)?;
            }
        }
    }

    Ok(())
}
