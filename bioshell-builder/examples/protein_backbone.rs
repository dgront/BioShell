use std::env;
use std::path::{Path};
use bioshell_builder::{BuilderError, InternalCoordinatesDatabase, KinematicAtomTree};

const N_RESIDUES: usize = 10;

fn main() -> Result<(), BuilderError>{
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();
    let args: Vec<String> = env::args().collect();

    // --- build the database of internal monomers' definitions
    let bioshell_cif_monomers = Path::new("../../../bioshell-builder/data/");
    let db_path = match Path::new(&args[0]).parent() {
        None => { Path::new("../../../bioshell-builder/data/").to_path_buf() }
        Some(path) => { path.join(bioshell_cif_monomers).canonicalize().unwrap() }
    };
    let mut db = InternalCoordinatesDatabase::from_cif_directory(db_path.to_str().unwrap())
        .expect(format!("Can't find a folder with CIF files to load: {db_path:?}").as_str());

    // --- create a chain builder
    let mut bb_builder = KinematicAtomTree::new();
    let bb_def = db.get_definition("bb_").unwrap();
    let cterm_def = db.get_definition("patch_CTerm").unwrap();

    // --- build a backbone chain
    for i in 0..N_RESIDUES {
        bb_builder.add_residue(&bb_def);
    }
    bb_builder.patch_residue(N_RESIDUES-1, &cterm_def)?;
    let atoms = bb_builder.build_atoms("A")?;

    // for iatom in 0..atoms.len() {
    //     let res_seq = bb_builder.residue_for_atom((iatom));
    //     let elem = bb_builder.atom_element(iatom).clone();
    //     println!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           {:2}",
    //              iatom + 1, &bb_builder.atom_name(iatom),
    //              "A", res_seq, &atoms[iatom].x, &atoms[iatom].y, &atoms[iatom].z, elem);
    // }

    return Ok(());
}
