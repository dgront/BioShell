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

/// Computes the Phi dihedral angle for a given amino acid residue.
///
/// Returns a [`PDBError::NoSuchAtom`](PDBError::NoSuchAtom) if the requested residue
/// or the preceding one misses any of the required atoms. Also returns
/// [`PDBError::ResidueAtChainBreak`](PDBError::ResidueAtTerminus) it the preceding residue can't be found.
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
/// let strctr = Structure::from_atoms("1xyz", atoms);
/// let phi_val = phi(&strctr, &ResidueId::new("A", 23, ' ')).unwrap();
/// assert_delta!(phi_val.to_degrees(), -71.925, 0.01);
/// ```
pub fn phi(strctr: &Structure, which_res: &ResidueId) -> Result<f64, PDBError> {
    let i_residue = strctr.residue_pos(which_res)?;
    if i_residue == 0 {
        return Err(PDBError::ResidueAtTerminus { res_id: which_res.clone() });
    }
    let prev_res = strctr.residue_ids[i_residue-1].clone();
    if prev_res.chain_id != which_res.chain_id {
        return Err(PDBError::ResidueAtTerminus { res_id: which_res.clone() });
    }
    let c_prev = strctr.atom(&prev_res," C  ")?;
    let n = strctr.atom(&which_res," N  ")?;
    let ca = strctr.atom(&which_res," CA ")?;
    let c = strctr.atom(&which_res," C  ")?;

    return Ok(dihedral_angle4(&c_prev.pos, &n.pos, &ca.pos, &c.pos));
}

/// Computes the Psi dihedral angle for a given amino acid residue.
///
/// Returns a [`PDBError::NoSuchAtom`](PDBError::NoSuchAtom) if the requested residue
/// or the following one misses any of the required atoms. Also returns
/// [`PDBError::ResidueAtChainBreak`](PDBError::ResidueAtTerminus) it the following residue can't be found.
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
/// let strctr = Structure::from_atoms("1xyz", atoms);
/// let psi_val = psi(&strctr, &ResidueId::new("A", 23, ' ')).unwrap();
/// assert_delta!(psi_val.to_degrees(), -31.158, 0.01);
/// ```
pub fn psi(strctr: &Structure, which_res: &ResidueId) -> Result<f64, PDBError> {
    let i_residue = strctr.residue_pos(which_res)?;
    if i_residue == strctr.atoms_for_residue_id.len() - 1 {
        return Err(PDBError::ResidueAtTerminus { res_id: which_res.clone() });
    }
    let next_res = strctr.residue_ids[i_residue + 1].clone();
    if next_res.chain_id != which_res.chain_id {
        return Err(PDBError::ResidueAtTerminus { res_id: which_res.clone() });
    }
    let n = strctr.atom(&which_res," N  ")?;
    let ca = strctr.atom(&which_res," CA ")?;
    let c = strctr.atom(&which_res," C  ")?;
    let n_next = strctr.atom(&next_res," N  ")?;

    return Ok(dihedral_angle4(&n.pos, &ca.pos, &c.pos, &n_next.pos));
}