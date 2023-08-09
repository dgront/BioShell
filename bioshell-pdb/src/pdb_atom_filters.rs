//! Allows to process atoms of a [`Structure`](crate::Structure) with iterators.
//!
//! Atoms of a [`Structure`](crate::Structure) may be accessed by [`Structure::atoms()`](crate::Structure::atoms())
//! or [`Structure::atoms_mut()`](crate::Structure::atoms_mut())
//! (immutable or mutable access, respectively) which borrow a reference to a vector of atoms.
//! These can be processed in a standard rust way with iterator tools. A [`PdbAtomPredicate`](PdbAtomPredicate)
//! trait defines [`check()`](PdbAtomPredicate::check()) method that returns true if a given predicate is satisfied.
//!
//! For example, [`IsBackbone`](IsBackbone) predicate may be used to filter backbone atoms:
//! ```
//! use bioshell_pdb::{PdbAtom, Structure};
//! use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! let mut strctr = Structure::new();
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! let bb = IsBackbone{};
//! for bb_atom in strctr.atoms().iter().filter(|b| bb.check(b)) {
//!     println!("{}",bb_atom.name);
//! }
//!```
//!
//! Atoms selected from a [`Structure`](crate::Structure) may be **cloned and collected** into a new `Vec<PdbAtom>` as below:
//! ```
//! # use bioshell_pdb::{PdbAtom, Structure};
//! # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! # let mut strctr = Structure::new();
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! let bb = IsBackbone{};
//! let bb_only: Vec<PdbAtom> = strctr.atoms().iter().filter(|b| bb.check(b)).cloned().collect();
//! assert_eq!(bb_only.len(), 4);
//! ```
//! Alternatively, an iteration over atoms may be directly copied into a new [`Structure`](crate::Structure):
//! ```
//! # use bioshell_pdb::{PdbAtom, Structure};
//! # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! # let mut strctr = Structure::new();
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! # let bb = IsBackbone{};
//! let bb_strctr = Structure::from_iterator(strctr.atoms().iter().filter(|b| bb.check(b)));
//! ```

use crate::{PdbAtom, ResidueId};

/// A handy filter to process atoms of a [`Structure`](crate::Structure) with iterators.
///
/// Structs implementing [`PdbAtomPredicate`](PdbAtomPredicate) trait can be used as predicates
/// while filtering Rust iterators. Example below shows how to iterate over backbone atoms
pub trait PdbAtomPredicate {
    fn check(&self, a: &PdbAtom) -> bool;
}

/// Returns `true` if an atom belongs to a certain chain.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByChain, PdbAtomPredicate};
/// let mut strctr = Structure::new();
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_chain_A = ByChain::new("A");
/// let atoms_A: Vec<PdbAtom> = strctr.atoms().iter().filter(|a| select_chain_A.check(a)).cloned().collect();
/// assert_eq!(atoms_A.len(), 1);
/// ```
pub struct ByChain {chain_id: String}

impl ByChain {
    pub fn new(code: &str) -> Self { ByChain { chain_id: String::from(code) } }
}

impl PdbAtomPredicate for ByChain {
    fn check(&self, a: &PdbAtom) -> bool { a.chain_id == self.chain_id }
}

/// Returns `true` if an atom belongs to a certain residue.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByChain, PdbAtomPredicate};
/// let mut strctr = Structure::new();
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_chain_A = ByChain::new("A");
/// let atoms_A: Vec<PdbAtom> = strctr.atoms().iter().filter(|a| select_chain_A.check(a)).cloned().collect();
/// assert_eq!(atoms_A.len(), 1);
/// ```
pub struct ByResidue {res_id: ResidueId }

impl ByResidue {
    pub fn new(res_id: ResidueId) -> ByResidue { ByResidue { res_id } }
}

impl PdbAtomPredicate for ByResidue {
    fn check(&self, a: &PdbAtom) -> bool { self.res_id.check(a) }
}

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
pub struct IsBackbone;

impl PdbAtomPredicate for IsBackbone {
    fn check(&self, a: &PdbAtom) -> bool {
        a.name == " CA " || a.name == " C  " || a.name == " N  " || a.name == " O  "
    }
}

/// Returns `true` if an atom is an alpha carbon
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
/// # let mut strctr = Structure::new();
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
/// let ca = IsCA{};
/// let ca_count = strctr.atoms().iter().filter(|a| ca.check(a)).count();
/// # assert_eq!(ca_count, 1);
/// ```
pub struct IsCA;

impl PdbAtomPredicate for IsCA {
    fn check(&self, a: &PdbAtom) -> bool { a.name == " CA " }
}

/// A filter defined for a pair of atoms.
///
/// Structs implementing [`PdbAtomPredicate`](PdbAtomPredicate) trait can be used as predicates
/// while filtering Rust iterators. Example below shows how to iterate over backbone atoms
pub trait PdbAtomPredicate2 {
    fn check(&self, ai: &PdbAtom, aj: &PdbAtom) -> bool;
}

/// Returns `true` if both atoms belong to the same chain
pub struct SameChain;

impl PdbAtomPredicate2 for SameChain {
    fn check(&self, ai: &PdbAtom, aj: &PdbAtom) -> bool { ai.chain_id == aj.chain_id }
}

/// Returns `true` if both atoms belong to the same residue
pub struct SameResidue;

impl PdbAtomPredicate2 for SameResidue {
    fn check(&self, ai: &PdbAtom, aj: &PdbAtom) -> bool {
        ai.chain_id == aj.chain_id && ai.res_seq == aj.res_seq && ai.i_code == aj.i_code
    }
}