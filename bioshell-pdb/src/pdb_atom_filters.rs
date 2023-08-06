use crate::PdbAtom;

/// A handy filter to process atoms of a  [`Structure`](Structure) with iterators.
///
/// # Examples
/// ```
/// use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
/// let mut strctr = Structure::new();
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
/// let bb = IsBackbone{};
/// let bb_count = strctr.atoms().iter().filter(|b| bb.check(b)).count();
/// assert_eq!(bb_count, 4);
/// ```
pub trait PdbAtomPredicate {
    fn check(&self, a: &PdbAtom) -> bool;
}

/// Returns `true` if an atom belongs to a certain chain.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ChainMatches, PdbAtomPredicate};
/// let mut strctr = Structure::new();
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_chain_A = ChainMatches::new("A");
/// let atoms_A: Vec<PdbAtom> = strctr.atoms().iter().filter(|a| select_chain_A.check(a)).cloned().collect();
/// assert_eq!(atoms_A.len(), 1);
/// ```
pub struct ChainMatches {chain_id: String}

impl ChainMatches {
    pub fn new(code: &str) -> Self { ChainMatches{ chain_id: String::from(code) } }
}

impl PdbAtomPredicate for ChainMatches {
    fn check(&self, a: &PdbAtom) -> bool { a.chain_id == self.chain_id }
}

pub struct IsBackbone;

/// Returns `true` if an atom belongs to a backbone.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
/// # let mut strctr = Structure::new();
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
/// let bb = IsBackbone{};
/// let bb_count = strctr.atoms().iter().filter(|b| bb.check(b)).count();
/// # assert_eq!(bb_count, 4);
/// ```
impl PdbAtomPredicate for IsBackbone {
    fn check(&self, a: &PdbAtom) -> bool {
        a.name == " CA " || a.name == " C  " || a.name == " N  " || a.name == " O  "
    }
}



