use std::env;
use bioshell_chem::{load_molecule};
use bioshell_chem::icoords::{ZMatrix};


fn main() -> Result<(), Box<dyn std::error::Error>> {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    // ---------- Prepare the molecule and the kinematic atom chain
    let mut molcle = load_molecule(&args[1])?;
    let zmatrix = ZMatrix::from_molecule(&mut molcle, 0, 1, 2)?;
    println!("{}", &zmatrix);

    Ok(())
}
