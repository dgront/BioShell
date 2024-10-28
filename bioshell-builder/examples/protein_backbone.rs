use std::env;
use std::path::{Path};
use std::str::FromStr;
use std::time::Instant;
use log::info;
use rand::Rng;
use bioshell_builder::{BuilderError, InternalCoordinatesDatabase, KinematicAtomTree};

const N_RESIDUES: usize = 100;
const N_MODELS: usize = 1000;
const CONFORMATION: Conformation = Conformation::Helical(helical_phi_psi);
enum Conformation {
    Random(fn() -> (f64, f64)),
    Helical(fn() -> (f64, f64)),
}

// Implementing FromStr for Conformation
impl FromStr for Conformation {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "random" => Ok(Conformation::Random(random_phi_psi)),
            "helical" => Ok(Conformation::Helical(helical_phi_psi)),
            _ => Err(format!("Invalid conformation type: '{}'", s)),
        }
    }
}

fn random_phi_psi() -> (f64, f64) {
    let mut rng = rand::thread_rng();
    (rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI), rng.gen_range(-std::f64::consts::PI..std::f64::consts::PI))
}

fn helical_phi_psi() -> (f64, f64) {
    let mut rng = rand::thread_rng();
    let sd = 0.36_f64.to_radians();
    (-64.0_f64.to_radians() + rng.gen_range(-sd..sd), -41.0_f64.to_radians() + rng.gen_range(-sd..sd))
}

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

    // --- modify Phi, Psi angles
    let phi_psi_gen = match CONFORMATION {
        Conformation::Random(f) => {f}
        Conformation::Helical(f) => {f}
    };

    let start = Instant::now();
    for i_model in 0..N_MODELS {
        println!("MODEL    {:3}", i_model + 1);
        for i in 0..N_RESIDUES {
            let (phi, psi) = phi_psi_gen();
            bb_builder.set_named_dihedral(i, "Phi", phi)?;
            bb_builder.set_named_dihedral(i, "Psi", psi)?;
        }
        let atoms = bb_builder.build_atoms("A")?;
        for mut a in atoms {
            a.res_name = "GLY".to_string();
            println!("{}", a);
        }
        println!("ENDMDL");
    }
    info!("{} chains built in: {:?}", N_MODELS, start.elapsed());

    return Ok(());
}
