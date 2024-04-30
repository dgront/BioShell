use crate::{PdbAtom, ResidueId, Structure};
use crate::calc::dihedral_angle4;
use crate::pdb_parsing_error::PDBError;

/// Calculate the squared distance between two atoms.
pub fn distance_squared(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    return ai.pos.distance_square_to(&aj.pos);
}

/// Calculate the distance between two atoms.
///
/// # Example
///```rust
/// use bioshell_pdb::PdbAtom;
/// use bioshell_pdb::calc::distance;
/// let ai = PdbAtom::from_atom_line("ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C");
/// let aj = PdbAtom::from_atom_line("ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C");
/// let d = distance(&ai, &aj);
/// # assert!((d-3.8).abs()<0.1);
///```
pub fn distance(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    return ai.pos.distance_to(&aj.pos);
}

///
/// # Example
/// ```
/// use bioshell_pdb::{assert_delta, PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::calc::phi;
/// let pdb_lines = ["ATOM    339  C   ASP A  22     -10.767   3.482  -3.375  1.00  0.20           C",
/// "ATOM    349  N   ALA A  23     -10.594   2.767  -4.455  1.00  0.19           N",
/// "ATOM    350  CA  ALA A  23      -9.237   2.226  -4.777  1.00  0.18           C",
/// "ATOM    351  C   ALA A  23      -8.287   3.342  -5.231  1.00  0.19           C"];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_iterator("1xyz", atoms.iter());
/// let phi_val = phi(&strctr, &ResidueId::new("A", 23, ' ')).unwrap();
/// assert_delta!(phi_val.to_degrees(), -71.925, 0.01);
/// ```
pub fn phi(strctr: &Structure, which_res: &ResidueId) -> Result<f64, PDBError> {
    let i_residue = strctr.residue_pos(which_res)?;
    if i_residue == 0 { return Ok(0.0); }
    let prev_res = strctr.residue_ids[i_residue-1].clone();
    let c_prev = strctr.atom(&prev_res," C  ")?;
    let n = strctr.atom(&which_res," N  ")?;
    let ca = strctr.atom(&which_res," CA ")?;
    let c = strctr.atom(&which_res," C  ")?;

    return Ok(dihedral_angle4(&c_prev.pos, &n.pos, &ca.pos, &c.pos));
}

/// Computes the Psi dihedral angle for a given amino acid residue.
///
/// # Example
/// ```
/// use bioshell_pdb::{assert_delta, PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::calc::psi;
/// let pdb_lines = ["ATOM    349  N   ALA A  23     -10.594   2.767  -4.455  1.00  0.19           N",
/// "ATOM    350  CA  ALA A  23      -9.237   2.226  -4.777  1.00  0.18           C",
/// "ATOM    351  C   ALA A  23      -8.287   3.342  -5.231  1.00  0.19           C",
/// "ATOM    359  N   ALA A  24      -8.828   4.348  -5.860  1.00  0.22           N "];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_iterator("1xyz", atoms.iter());
/// let psi_val = psi(&strctr, &ResidueId::new("A", 23, ' ')).unwrap();
/// assert_delta!(psi_val.to_degrees(), -31.158, 0.01);
/// ```
pub fn psi(strctr: &Structure, which_res: &ResidueId) -> Result<f64, PDBError> {
    let i_residue = strctr.residue_pos(which_res)?;
    if i_residue == strctr.atoms_for_residueid.len() - 1 { return Ok(0.0); }
    let next_res = strctr.residue_ids[i_residue + 1].clone();
    let n = strctr.atom(&which_res," N  ")?;
    let ca = strctr.atom(&which_res," CA ")?;
    let c = strctr.atom(&which_res," C  ")?;
    let n_next = strctr.atom(&next_res," N  ")?;

    return Ok(dihedral_angle4(&n.pos, &ca.pos, &c.pos, &n_next.pos));
}