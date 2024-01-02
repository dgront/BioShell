use std::env;
use std::io::BufReader;
use bioshell_builder::{BuilderError, InternalCoordinatesDatabase, KinematicAtomTree};
use bioshell_cif::read_cif_buffer;

const N_RESIDUES: usize = 10;

fn main() -> Result<(), BuilderError>{
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    // --- build the database of internal monomers' definitions
    let mut db = InternalCoordinatesDatabase::new();
    let cif_data = read_cif_buffer(BufReader::new(BB_.as_bytes()));
    db.load_from_cif_data(cif_data);
    let cif_data = read_cif_buffer(BufReader::new(CTerm.as_bytes()));
    db.load_from_cif_data(cif_data);

    // --- create a chain builder
    let mut bb_builder = KinematicAtomTree::new();
    let bb_def = db.get_definition("bb_").unwrap();
    let cterm_def = db.get_definition("patch_CTerm").unwrap();

    // --- build a backbone chain
    for i in 0..N_RESIDUES {
        bb_builder.add_residue(&bb_def);
    }
    bb_builder.patch_residue(N_RESIDUES-1, &cterm_def);
    let atoms = bb_builder.restore_atoms()?;

    for iatom in 0..atoms.len() {
        let res_seq = bb_builder.residue_for_atom((iatom));
        println!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C",
                 iatom + 1, &bb_builder.atom_name(iatom),
                 "A", res_seq, &atoms[iatom].x, &atoms[iatom].y, &atoms[iatom].z);
    }

    return Ok(());
}


const BB_: &str = "data_bb_
loop_
_res_name
_atom_a_residue_locator
_atom_a_name
_atom_b_residue_locator
_atom_b_name
_atom_c_residue_locator
_atom_c_name
_atom_d_residue_locator
_atom_d_name
_c_d_bond_length
_b_c_d_planar_angle
_a_b_c_d_dihedral_angle
_dihedral_angle_name
'bb ' prev ' N  ' prev ' CA ' prev ' C  ' this ' N  ' 1.328685 114.0  180.0 psi
'bb ' prev ' CA ' prev ' C  ' this ' N  ' this ' CA ' 1.458001 123.0  180.0 omega
'bb ' prev ' C  ' this ' N  ' this ' CA ' this ' C  ' 1.523258 110.0 -180.0 phi
'bb ' next ' N  ' this ' CA ' this ' C  ' this ' O  ' 1.231015 120.0  180.0 -
#
";

const CTerm: &str = "data_patch_CTerm
loop_
_res_name
_atom_a_residue_locator
_atom_a_name
_atom_b_residue_locator
_atom_b_name
_atom_c_residue_locator
_atom_c_name
_atom_d_residue_locator
_atom_d_name
_c_d_bond_length
_b_c_d_planar_angle
_a_b_c_d_dihedral_angle
_dihedral_angle_name
'CTerm' this ' N  ' this ' CA ' this ' C  ' this ' OXT' 1.2      116.5  180.0 psi
'CTerm' this ' OXT' this ' CA ' this ' C  ' this ' O  ' 1.231015 120.0  180.0 -
#
";