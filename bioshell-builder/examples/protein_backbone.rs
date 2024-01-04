use std::env;
use std::path::{Path};
use rand::Rng;
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
    let db = InternalCoordinatesDatabase::from_cif_directory(db_path.to_str().unwrap())
        .expect(format!("Can't find a folder with CIF files to load: {db_path:?}").as_str());

    // --- create a chain builder
    let mut bb_builder = KinematicAtomTree::new();
    let bb_def = db.get_definition("bb_").unwrap();
    let cterm_def = db.get_definition("patch_CTerm").unwrap();

    // --- build a backbone chain
    for _i in 0..N_RESIDUES {
        bb_builder.add_residue(&bb_def);
    }
    bb_builder.patch_residue(N_RESIDUES-1, &cterm_def)?;
    // --- randomize Phi, Psi angles
    let mut rng = rand::thread_rng();
    for i in 0..N_RESIDUES {
        let phi: f64 = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
        let psi: f64 = rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI);
        bb_builder.set_named_dihedral(i, "Phi", phi)?;
        bb_builder.set_named_dihedral(i, "Psi", psi)?;
    }

    let atoms = bb_builder.build_atoms("A")?;
    for a in atoms { println!("{}", a); }

    return Ok(());
}
