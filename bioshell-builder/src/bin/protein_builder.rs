use std::env;
use clap::{Parser};
use bioshell_builder::{BuilderError, InternalCoordinatesDatabase, KinematicAtomTree};
use bioshell_cif::read_cif_buffer;
use bioshell_io::open_file;
use bioshell_pdb::{PdbAtom};
use bioshell_seq::sequence::{FastaIterator, load_sequences};
use bioshell_seq::SequenceError;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Command line tool to operate on PDB files
/// say pdb_tool -h to see options
struct Args {
    /// output PDB file name
    #[clap(short, long, short='o', required=true)]
    outfile: String,
    /// input protein sequence in the FASTA format
    #[clap(short, long, short='f')]
    infasta: String,
}

fn build(aa_codes: Vec<String>, chain_id: &str) -> Result<Vec<PdbAtom>, BuilderError> {
    let mut output: Vec<PdbAtom> = vec![];

    let mut residue_manager = InternalCoordinatesDatabase::new();
    let mut cif_reader = open_file("./data/bb_topo.cif")?;
    let data_blocks = read_cif_buffer(&mut cif_reader);
    residue_manager.load_from_cif_data(data_blocks);

    // --- create a chain builder
    let mut bb_builder = KinematicAtomTree::new();
    let bb_def = residue_manager.get_definition("bb_").unwrap();
    let cterm_def = residue_manager.get_definition("patch_CTerm").unwrap();
    bb_builder.add_residue(&bb_def);
    bb_builder.add_residue(&bb_def);
    bb_builder.patch_residue(1, &cterm_def)?;
    let atoms = bb_builder.build_atoms("A")?;

    return Ok(atoms);
}


fn main() -> Result<(), SequenceError> {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // ---------- INPUT section
    let seq = load_sequences(&args.infasta, "target")?;
    let mut sequences: Vec<(String, String)> = vec![];
    for sequence in &seq {
        let id = String::from(sequence.id());
        sequences.push((sequence.to_string(), id));
    }

    Ok(())
}

