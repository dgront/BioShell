//! Provides filters to process atoms of a [`Structure`](crate::Structure) with iterators.
//!
//! Each predicate (i.e. a filter) provided by this module, implements [`PdbAtomPredicate`](PdbAtomPredicate) trait.
//! The [`PdbAtomPredicate::check()`](PdbAtomPredicate::check()) method returns true if a given predicate is satisfied.
//! This facilitates filtering [`PdbAtom`](crate::PdbAtom)s of a given [`Structure`](crate::Structure)
//! in a standard rust way with iterator tools.
//! Atoms of a [`Structure`](crate::Structure) may be accessed by [`Structure::atoms()`](crate::Structure::atoms())
//! or [`Structure::atoms_mut()`](crate::Structure::atoms_mut())
//! (immutable or mutable access, respectively), which borrow a reference to a vector of atoms.
//!
//! For example, [`IsBackbone`](IsBackbone) predicate may be used to filter backbone atoms:
//! ```
//! use bioshell_pdb::{PdbAtom, Structure};
//! use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! let mut strctr = Structure::new("1xyz");
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! let bb = IsBackbone;
//! for bb_atom in strctr.atoms().iter().filter(|b| bb.check(b)) {
//!     println!("{}",bb_atom.name);
//! }
//!```
//!
//! Atoms selected from a [`Structure`](crate::Structure) may be **cloned and collected** into a new `Vec<PdbAtom>` as below:
//! ```
//! # use bioshell_pdb::{PdbAtom, Structure};
//! # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! # let mut strctr = Structure::new("1xyz");
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! let bb = IsBackbone;
//! let bb_only: Vec<PdbAtom> = strctr.atoms().iter().filter(|b| bb.check(b)).cloned().collect();
//! assert_eq!(bb_only.len(), 4);
//! ```
//! Alternatively, an iteration over atoms may be directly copied into a new [`Structure`](crate::Structure):
//! ```
//! # use bioshell_pdb::{PdbAtom, Structure};
//! # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
//! # let mut strctr = Structure::new("1xyz");
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
//! # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
//! # let bb = IsBackbone;
//! let bb_strctr = Structure::from_iterator("1xyz", strctr.atoms().iter().filter(|b| bb.check(b)));
//! ```

use bioshell_seq::chemical::ResidueType;
use crate::{PdbAtom, ResidueId};

/// A handy filter to process atoms of a [`Structure`](crate::Structure) with iterators.
///
/// Structs implementing [`PdbAtomPredicate`](PdbAtomPredicate) trait can be used as predicates
/// while filtering Rust iterators. Example below shows how to iterate over backbone atoms
pub trait PdbAtomPredicate {
    fn check(&self, a: &PdbAtom) -> bool;
}

/// Always returns `true`
///
/// Declare this [`PdbAtomPredicate`](PdbAtomPredicate) when you have to use one but you don't want
/// to skip any atoms
pub struct AlwaysPass;

impl PdbAtomPredicate for AlwaysPass {
    fn check(&self, _a: &PdbAtom) -> bool {true}
}

/// Returns `true` if an atom belongs to a certain chain.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByChain, PdbAtomPredicate};
/// let mut strctr = Structure::new("1xyz");
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
/// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByResidue, PdbAtomPredicate};
/// let mut strctr = Structure::new("1xyz");
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_ala_A = ByResidue::new(ResidueId::new("A", 69, ' '));
/// let ala_A: Vec<PdbAtom> = strctr.atoms().iter().filter(|a| select_ala_A.check(a)).cloned().collect();
/// assert_eq!(ala_A.len(), 1);
/// ```
pub struct ByResidue {res_id: ResidueId }

impl ByResidue {
    pub fn new(res_id: ResidueId) -> ByResidue { ByResidue { res_id } }
}

/// Returns `true` if an atom belongs to a certain residue.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByResidue, PdbAtomPredicate};
/// let mut strctr = Structure::new("1xyz");
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_ala_A = ByResidue::new(ResidueId::new("A", 69, ' '));
/// let ala_A: Vec<PdbAtom> = strctr.atoms().iter().filter(|a| select_ala_A.check(a)).cloned().collect();
/// assert_eq!(ala_A.len(), 1);
/// ```
pub struct ByResidueType {res_type: ResidueType }

impl ByResidue {
    pub fn new(res_id: ResidueId) -> ByResidue { ByResidue { res_id } }
}

impl PdbAtomPredicate for ByResidue {
    fn check(&self, a: &PdbAtom) -> bool { self.res_id.check(a) }
}

/// Returns `true` if an atom belongs to a given residue range.
///
/// The range is defined by two [`ResidueId`](crate::ResidueId)s objects, both inclusive.
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
/// use bioshell_pdb::pdb_atom_filters::{ByResidueRange, PdbAtomPredicate};
/// let mut strctr = Structure::new("1xyz");
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  67      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  CA  ALA A  68A     25.790  28.757  29.513  1.00 16.12           C"));
/// strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// let select_ala_A = ByResidueRange::new(ResidueId::new("A", 68, ' '),ResidueId::new("A", 68, 'A'));
/// let cnt = strctr.atoms().iter().filter(|a| select_ala_A.check(a)).count();
/// assert_eq!(cnt, 2);
/// ```
pub struct ByResidueRange { first_res_id: ResidueId, last_res_id: ResidueId }

impl ByResidueRange {
    pub fn new(first_res_id: ResidueId, last_res_id: ResidueId ) -> ByResidueRange {
        ByResidueRange { first_res_id, last_res_id }
    }
}

impl PdbAtomPredicate for ByResidueRange {
    fn check(&self, a: &PdbAtom) -> bool {
        if a.chain_id < self.first_res_id.chain_id || a.chain_id > self.last_res_id.chain_id { return false }
        if (a.res_seq < self.first_res_id.res_seq && a.chain_id==self.first_res_id.chain_id)
            || (a.res_seq > self.last_res_id.res_seq && a.chain_id==self.last_res_id.chain_id) { return false }
        if a.i_code == ' ' {
            if self.first_res_id.i_code != ' '
                && self.first_res_id.chain_id == a.chain_id
                && self.first_res_id.res_seq == a.res_seq { return false }
        } else {
            if self.last_res_id.chain_id == a.chain_id && self.last_res_id.res_seq == a.res_seq {
                if a.i_code < self.first_res_id.i_code { return false }
                if self.last_res_id.i_code == ' ' { return false }
            }
        }

        return true;
    }
}

/// Returns `true` if an atom belongs to a backbone.
///
/// The predicate returns true for protein backbone heavy atoms: `N`, `CA`, `C`, `O`, `OXT` as well as
/// for hydrogens: `H`, `HA`, `HA2` and `HA3`.
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
/// # let mut strctr = Structure::new("1xyz");
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
/// let bb = IsBackbone;
/// let bb_count = strctr.atoms().iter().filter(|b| bb.check(b)).count();
/// # assert_eq!(bb_count, 4);
/// ```
pub struct IsBackbone;

impl PdbAtomPredicate for IsBackbone {
    fn check(&self, a: &PdbAtom) -> bool {
        a.name == " CA " || a.name == " C  " || a.name == " N  " || a.name == " O  " || a.name == " H  " || a.name == " OXT"
            || a.name == " HA " || a.name == " HA2" || a.name == " HA3"
    }
}

/// Returns `true` if an atom is an alpha carbon
///
/// # Examples
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
/// # let mut strctr = Structure::new("1xyz");
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
/// let ca = IsCA;
/// let ca_count = strctr.atoms().iter().filter(|a| ca.check(a)).count();
/// # assert_eq!(ca_count, 1);
/// ```
pub struct IsCA;

impl PdbAtomPredicate for IsCA {
    fn check(&self, a: &PdbAtom) -> bool { a.name == " CA " }
}

/// Returns `true` if an atom is a hydrogen.
///
/// The following example removes all hydrogen atoms from a structure
///
/// # Examples
/// ```
/// use bioshell_pdb::{load_pdb_reader, PdbAtom, Structure};
/// use std::io::BufReader;
/// use bioshell_pdb::pdb_atom_filters::{IsHydrogen, PdbAtomPredicate};
/// let gly_pdb = "ATOM    148  N   GLY A   9       9.692  -3.742   0.370  1.00  0.14           N
/// ATOM    149  CA  GLY A   9      10.920  -2.963   0.070  1.00  0.18           C
/// ATOM    150  C   GLY A   9      12.144  -3.871   0.156  1.00  0.20           C
/// ATOM    151  O   GLY A   9      12.210  -4.757   0.986  1.00  0.24           O
/// ATOM    152  H   GLY A   9       9.093  -3.471   1.097  1.00  0.14           H
/// ATOM    153  HA2 GLY A   9      10.848  -2.565  -0.927  1.00  0.20           H
/// ATOM    154  HA3 GLY A   9      11.018  -2.149   0.770  1.00  0.21           H";
/// let mut strctr = load_pdb_reader(BufReader::new(gly_pdb.as_bytes())).unwrap();
/// let is_h = IsHydrogen;
/// strctr.atoms_mut().retain(|a| !is_h.check(a));
/// # assert_eq!(strctr.count_atoms(), 4);
/// ```
pub struct IsHydrogen;

impl PdbAtomPredicate for IsHydrogen {
    fn check(&self, a: &PdbAtom) -> bool {
        if let Some(element) = &a.element {
            if element == "H" { return true }
        }
        return false;
    }
}

/// Returns `true` if an atom belongs to a water molecule
///
/// # Examples
/// The following example removes water molecules from a structure:
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsWater, PdbAtomPredicate};
/// # let mut strctr = Structure::new("1xyz");
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  O   HOH A  69      25.155  27.554  29.987  1.00 21.91           O"));
/// let hoh = IsWater;
/// strctr.atoms_mut().retain(|a| !hoh.check(&a));
/// # assert_eq!(strctr.count_atoms(), 1);
/// # assert_eq!(strctr.atoms()[0].name, " CA ");
/// ```
pub struct IsWater;

impl PdbAtomPredicate for IsWater {
    fn check(&self, a: &PdbAtom) -> bool { a.res_name == "HOH" }
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