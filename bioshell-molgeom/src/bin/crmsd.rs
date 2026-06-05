use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use clap::Parser;
use anyhow::{Context, Result};

use bioshell_pdb::{Deposit, PdbAtom};
use bioshell_pdb::pdb_atom_filters::{ByChain, InvertPredicate, IsBackbone, IsCA, IsHydrogen, KeepProtein, MatchAll, PdbAtomPredicate};
use log::info;
use bioshell_core::Vec3;
use bioshell_molgeom::align::crmsd_transform;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None, arg_required_else_help = true)]
/// Calculates crmsd distance between two protein structures.
/// say distance_map -h to see options
struct Args {
    /// the first input protein structure in either PDB or mmCIF format
    #[clap(short, long, short='a', required=true)]
    protein_a: String,
    /// the second input protein structure in either PDB or mmCIF format
    #[clap(short, long, short='b', required=true)]
    protein_b: String,
    /// superimpose structures by matching atom names; if not set, atoms will be matched by their order in the structure file
    #[clap(long)]
    match_names: bool,
    /// use only C-alpha atoms for calculations
    #[clap(long)]
    ca_only: bool,
    /// use only backbone atoms for calculations
    #[clap(long)]
    bb_only: bool,
    /// keep only selected chain; the same chain will be selected from both structures
    #[clap(long)]
    select_chain: Option<String>,
    #[clap(long, short='o')]
    /// save transformed structure A in the PDB format
    output: Option<String>,
    /// be more verbose and log program actions on the screen
    #[clap(long, short='v')]
    verbose: bool
}

/// Returns atom hash used to pair identical atoms from two structures. The hash is based on the atom name, residue number and chain ID.
fn atom_hash(atom: &PdbAtom) -> String {
    return format!("{}{}{}{}", atom.name, atom.res_seq, atom.i_code, atom.chain_id);
}

fn main() -> Result<()> {
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
    if args.bb_only {
        multi_filter.add_predicate(Box::new(IsBackbone));
        multi_filter.add_predicate(Box::new(InvertPredicate::new(IsHydrogen)));
    }
    if args.ca_only { multi_filter.add_predicate(Box::new(IsCA)); }
    // --- if no filters were added, add default ones to keep only protein heavy atoms
    if multi_filter.count_filters() == 0 {
        multi_filter.add_predicate(Box::new(KeepProtein));
        multi_filter.add_predicate(Box::new(InvertPredicate::new(IsHydrogen)));
        info!("Selecting all non-hydrogen atoms from all protein chains");
    }
    // --- if chain selection is requested, add the corresponding filter
    if let Some(chain_id) = args.select_chain {
        info!("Selecting only chain {}", &chain_id);
        multi_filter.add_predicate(Box::new(ByChain::new(&chain_id)));
    }

    let deposit_a = Deposit::from_file(&args.protein_a).unwrap();
    let strctr_a = deposit_a.structure().unwrap();
    let atoms_a = strctr_a.atoms().iter()
        .filter(|a| multi_filter.check(a)).collect::<Vec<_>>();

    let deposit_b = Deposit::from_file(&args.protein_b).unwrap();
    let strctr_b = deposit_b.structure().unwrap();
    let atoms_b = strctr_b.atoms().iter()
        .filter(|a| multi_filter.check(a)).collect::<Vec<_>>();

    let mut pos_a: Vec<Vec3> = vec![];
    let mut pos_b: Vec<Vec3> = vec![];
    if args.match_names {
        info!("Matching atoms by their names");
        let mut name_map_a = std::collections::HashMap::new();
        for a in &atoms_a {
            let hash = atom_hash(a);
            name_map_a.insert(hash, a.pos.clone());
        }

        for atom_b in &atoms_b {
            let hash = atom_hash(atom_b);
            if let Some(pos_a_i) = name_map_a.get(&hash) {
                pos_a.push(*pos_a_i);
                pos_b.push(atom_b.pos);
            }
        }
    } else {
        pos_a = atoms_a.iter().map(|a| a.pos.clone()).collect::<Vec<_>>();
        pos_b = atoms_b.iter().map(|a| a.pos.clone()).collect::<Vec<_>>();
    }

    let (rms, rt) = crmsd_transform(&pos_a, &pos_b);
    println!("{:6.3} on {} atoms", rms, pos_a.len());

    if let Some(output_file) = args.output {
        info!("Saving transformed structure A to {}", &output_file);
        let file = File::create(&output_file)
            .with_context(|| format!("failed to create output file: {}", &output_file))?;
        let mut writer = BufWriter::new(file);
        for atom_a in &atoms_a {
            let pos = rt.apply(&atom_a.pos);
            writeln!(writer, "ATOM  {:5} {:>4} {:>3} {}{:4}    {:8.3}{:8.3}{:8.3}",
                atom_a.serial, atom_a.name, atom_a.res_name, atom_a.chain_id, atom_a.res_seq, pos.x, pos.y, pos.z)
                .with_context(|| format!("failed to write to output file: {}", &output_file))?;
        }
    }

    return Ok(());
}
