use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::io::Write;
use std::ops::Range;

use itertools::{Itertools};
use bioshell_seq::chemical::{ResidueType, ResidueTypeManager, ResidueTypeProperties, KNOWN_RESIDUE_TYPES};
use bioshell_seq::sequence::Sequence;

use crate::pdb_atom::{PdbAtom, same_residue_atoms};
use crate::pdb_atom_filters::{SameResidue, PdbAtomPredicate, PdbAtomPredicate2, ByResidueRange};
use crate::pdb_parsing_error::PDBError;
use crate::pdb_parsing_error::PDBError::{NoSuchAtom, NoSuchResidue};
use crate::{ResidueId, SecondaryStructureTypes};
use crate::calc::Vec3;
use crate::PDBError::WrongAtomsNumberInModel;
use crate::secondary_structure::SecondaryStructure;


/// A biomacromolecular structure composed of [`PdbAtom`](PdbAtom) objects.
///
/// A [`Structure`](Structure) struct holds all atoms in a `Vec<PdbAtom>` container; its implementation
/// provides methods to look at them in various ways.
///
/// # Creating a [`Structure`](Structure)
/// Typically one gets a [`Structure`](Structure) object by loading it from a file in the PDB format:
/// ```no_run
/// # use bioshell_pdb::PDBError;
/// # fn main() -> Result<(), PDBError> {
/// use bioshell_pdb::{Deposit, PDBError, Structure};
/// let deposit = Deposit::from_file("2gb1.pdb")?;
/// let strctr = deposit.structure();
/// # Ok(())
/// # }
/// ```
/// A [`Structure`](Structure) can be also created from an [`Iterator`](Iterator) over [`PdbAtom`](PdbAtom)s:
/// ```
/// use bioshell_pdb::{PdbAtom, Structure};
/// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// let strctr = Structure::from_atoms("1xyz", atoms);
/// # assert_eq!(strctr.count_atoms(), 4);
/// ```
///
/// # Accessing its atoms
/// A [`Structure`](Structure) implements two methods that provide mutable  and immutable borrow
/// of the vector of its atoms: [`atoms()`](Structure::atoms()) and [`atoms_mut()`](Structure::atoms_mut()) respectively.
/// These can be easily filtered by [`PdbAtomPredicate`](PdbAtomPredicate) predicates provided
/// by [`pdb_atom_filters`](crate::pdb_atom_filters) module.
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
/// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
/// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
/// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
/// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
/// # let strctr = Structure::from_atoms("1xyz", atoms);
/// let is_ca = IsCA;
/// # let mut n_ca = 0;
/// for ca in  strctr.atoms().iter().filter(|a| is_ca.check(&a)) {
///     // --- process each alpha carbon here
/// # n_ca += 1;
/// }
/// # assert_eq!(n_ca, 2);
/// ```
/// # Removing atoms, residues or chains
/// This can be easily done using the [`retain()`](std::vec::Vec::retain()) method of a [`Vec`](std::vec::Vec) struct combined with a respective
/// [`PdbAtomPredicate`](PdbAtomPredicate). An example below **removes water molecules**:
///
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsNotWater, PdbAtomPredicate};
/// # let mut strctr = Structure::new("1xyz");
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  O   HOH A  69      25.155  27.554  29.987  1.00 21.91           O"));
/// let hoh = IsNotWater;
/// let new_strctr = Structure::from_iterator("1xyz", strctr.atoms().iter().filter(|a| !hoh.check(&a)).cloned());
/// # assert_eq!(new_strctr.count_atoms(), 1);
/// ```
/// Another example removes all hydrogen atoms:
/// ```
/// # use bioshell_pdb::{PdbAtom, Structure};
/// use bioshell_pdb::pdb_atom_filters::{IsHydrogen, PdbAtomPredicate};
/// # let mut strctr = Structure::new("1xyz");
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    149  CA  GLY A   9      10.920  -2.963   0.070  1.00  0.18           C"));
/// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    153  HA2 GLY A   9      10.848  -2.565  -0.927  1.00  0.20           H"));
/// let is_h = IsHydrogen;
/// let new_strctr = Structure::from_iterator("1xyz", strctr.atoms().iter().filter(|a| !is_h.check(&a)).cloned());
/// # assert_eq!(new_strctr.count_atoms(), 1);
/// ```
///
#[derive(Clone)]
pub struct Structure {
    /// Four-character PDB code of this deposit, such as `2GB1` or `4HHB`
    pub id_code: String,
    pub(crate) atoms: Vec<PdbAtom>,
    pub(crate) model_coordinates: Vec<Vec<Vec3>>,
    /// id of each residue, in the order of their appearance in the structure;
    pub(crate) residue_ids: Vec<ResidueId>,
    /// range of atoms that belong to i-th residue; order is the same as in `residue_ids`
    pub(crate) atoms_for_residue_id: Vec<Range<usize>>,
}

impl Structure {
    /// Create a new empty [`Structure`] that contains no atoms.
    pub fn new(id_code: &str) -> Self {
        Self {
            id_code: id_code.to_string(),
            atoms: vec![],
            model_coordinates: vec![],
            residue_ids: vec![],
            atoms_for_residue_id: vec![],
        }
    }

    /// Creates a new [`Structure`](Structure) by filling it with provided atoms
    ///
    /// Atoms will be consumed in the process.
    ///
    /// # Example
    ///```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
    /// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    ///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    ///let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let strctr = Structure::from_atoms("1xyz", atoms);
    /// assert_eq!(strctr.count_atoms(), 4);
    /// ```
    pub fn from_atoms(id_code: &str, atoms: Vec<PdbAtom>) -> Structure {

        let mut strctr = Structure::new(id_code);
        strctr.atoms = atoms;
        strctr.update();

        return strctr;
    }

    /// Creates a new [`Structure`](Structure) by filling it with atoms from an iterator.
    ///
    /// Typically, the atoms come from another [`Structure`](Structure) which you want to modify.
    /// The atoms will be consumed in the process. You may clone them beforehand as in the example below.
    ///
    /// # Example
    ///```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// # use bioshell_pdb::pdb_atom_filters::{IsBackbone, PdbAtomPredicate};
    /// let mut strctr = Structure::new("1xyz");
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
    /// let bb = IsBackbone{};
    /// let bb_strctr = Structure::from_iterator(&strctr.id_code, strctr.atoms().iter().filter(|a|bb.check(a)).cloned());
    /// # assert_eq!(bb_strctr.count_atoms(), 4);
    /// ```
    pub fn from_iterator<T>(id_code: &str, iter: T) -> Structure
    where T: IntoIterator<Item = PdbAtom> {

        let mut strctr = Structure::new(id_code);
        for a in iter { strctr.atoms.push(a) }
        strctr.update();

        return strctr;
    }

    /// Renumbers atoms and residues of a copied [`Structure`](Structure).
    pub fn renumbered_structure(&self) -> Structure {
        let mut strctr = Structure::new(&self.id_code);
        let mut atom_id = 0;
        let mut resid_id = 1;
        let mut last_resid = self.atoms[0].res_seq;
        let mut last_icode = self.atoms[0].i_code;
        let mut last_chain = self.atoms[0].chain_id.clone();
        for a in self.atoms.iter() {
            atom_id += 1;
            if a.res_seq != last_resid || a.i_code != last_icode {
                resid_id += 1;
                last_resid = a.res_seq;
                last_icode = a.i_code;
            }
            if a.chain_id != last_chain {
                last_chain = a.chain_id.clone();
                resid_id = 1;
            }
            let mut a = a.clone();
            a.res_seq = resid_id;
            a.serial = atom_id;
            a.i_code = ' ';
            strctr.atoms.push(a);
        }
        strctr.update();
        return strctr;
    }

    /// Pushes a given [`PdbAtom`](PdbAtom) at the end of this [`Structure`](Structure)
    ///
    /// The given atom object will be consumed in the process. This method is inefficient as after
    /// every call the internal structure of residue / chain indexes must be reconstructed.
    /// Use [`Structure::from_iterator()`](Structure::from_iterator()) to construct a [`Structure`](Structure)
    /// from a large number of atoms rather than try to insert the one-by-one.
    pub fn push_atom(&mut self, a: PdbAtom) {

        self.atoms.push(a);
        self.update();
    }

    /// Counts atoms of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new("1xyz");
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// assert_eq!(strctr.count_atoms(), 2);
    /// ```
    pub fn count_atoms(&self) -> usize { self.atoms.len() }

    /// Counts residues of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new("1xyz");
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
    ///
    /// assert_eq!(strctr.count_residues(), 2);
    /// ```
    pub fn count_residues(&self) -> usize {
        let same_res = SameResidue{};
        return self.atoms().windows(2).filter(|a| !same_res.check(&a[0], &a[1])).count() + 1
    }

    /// Counts chains of this [`Structure`](Structure)
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let mut strctr = Structure::new("1xyz");
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"));
    /// strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA B  68      25.790  28.757  29.513  1.00 16.12           C"));
    ///
    /// assert_eq!(strctr.count_chains(), 2);
    /// ```
    pub fn count_chains(&self) -> usize {
        self.atoms.iter().map(|a| &a.chain_id).collect::<HashSet<_>>().len()
    }

    /// Counts models i.e. distinct conformations of this [`Structure`](Structure)
    pub fn count_models(&self) -> usize { self.model_coordinates.len() }

    pub fn set_model(&mut self, i_model: usize) -> Result<(), PDBError> {
        if self.model_coordinates[i_model].len() != self.atoms.len() {
            return Err(WrongAtomsNumberInModel { model_index: i_model });
        }
        for i in 0..self.model_coordinates[i_model].len() {
            self.atoms[i].pos.set(&self.model_coordinates[i_model][i]);
        }
        return Ok(());
    }

    /// Provides immutable access to an atom
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure, ResidueId};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_atoms("1xyz", atoms);
    /// let a = strctr.atom(&ResidueId::new("A", 69, ' ')," CA ").unwrap();
    /// assert_eq!(a.name, " CA ");
    /// assert_eq!(a.res_seq, 69);
    /// ```
    pub fn atom(&self, res_id: &ResidueId, name: &str) -> Result<&PdbAtom, PDBError> {
        let i_residue = self.residue_pos(res_id)?;
        let range = &self.atoms_for_residue_id[i_residue];
        for i in range.start..range.end {
            if self.atoms[i].name == name { return Ok(&self.atoms[i]) }
        }
        return Err(NoSuchAtom { atom_name: name.to_string(), res_id: res_id.clone() });
    }

    /// Provides mutable access to an atom
    pub fn atom_mut(&mut self, res_id: &ResidueId, name: &str) -> Result<&mut PdbAtom, PDBError> {
        let i_residue = self.residue_pos(res_id)?;
        let range = &self.atoms_for_residue_id[i_residue];
        for i in range.start..range.end {
            if self.atoms[i].name == name { return Ok(&mut self.atoms[i]) }
        }
        return Err(NoSuchAtom { atom_name: name.to_string(), res_id: res_id.clone() });
    }

    /// Provides immutable access to atoms of this [`Structure`](Structure)
    pub fn atoms(&self) -> &Vec<PdbAtom> { &self.atoms }

    /// Creates a vector that holds string identifiers for all chains of this [`Structure`](Structure)
    ///
    /// The vector is sorted alphabetically, regardless the order of chains in this [`Structure`](Structure)
    pub fn chain_ids(&self) -> Vec<String> {
        let uniq: HashSet<&String> = self.atoms.iter().map(|a| &a.chain_id).collect();
        let mut output = Vec::from_iter(uniq.iter().map(|s| *s).cloned());
        output.sort();

        return output;
    }

    /// Creates a vector that holds string identifiers for all entities of this [`Structure`](Structure)
    ///
    /// The vector is sorted alphabetically.
    pub fn entity_ids(&self) -> Vec<String> {
        let uniq: HashSet<&String> = self.atoms.iter().map(|a| &a.entity_id).collect();
        let mut output = Vec::from_iter(uniq.iter().map(|s| *s).cloned());
        output.sort();

        return output;
    }

    /// Returns atoms of a given chain
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_atoms("1xyz", atoms);
    /// # #[allow(non_snake_case)]
    /// let chain_A_atoms = strctr.atoms_in_chain("A");
    /// # assert_eq!(chain_A_atoms.count(),4);
    /// ```
    pub fn atoms_in_chain<'a>(&'a self, chain_id: &'a str) -> impl Iterator<Item = &'a PdbAtom> + 'a {
        self.atoms.iter().filter(move |atm| atm.chain_id == chain_id)
    }

    /// Returns atoms of a given residue
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// # let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    /// #                     "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    /// #                     "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// # let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// # let strctr = Structure::from_atoms("1xyz", atoms);
    /// let res_atoms = strctr.atoms_in_residue(&ResidueId::new("A", 68, ' ')).unwrap();
    /// # assert_eq!(res_atoms.count(),2);
    /// ```
    pub fn atoms_in_residue(&self, residue_id: &ResidueId) -> Result<impl Iterator<Item = &PdbAtom>, PDBError> {
        if let Some(pos) = self.residue_ids.iter().position(|x| x == residue_id) {
            return Ok(self.atoms_for_residue_id[pos].clone().map(|i| &self.atoms[i]));
        } else {
            return  Err(NoSuchResidue{res_id: residue_id.clone()})
        }
    }

    /// Provide an iterator over atoms from a given range of residues
    ///
    /// The atoms may belong to different chains.
    ///
    /// # Example
    /// ```
    /// # use std::io::BufReader;
    /// use bioshell_pdb::{Deposit, ResidueId};
    /// # #[allow(non_upper_case_globals)]
    /// const pdb_txt: &str =
    /// "ATOM      2  CA  MET A   1     -13.296   0.028   3.924  1.00  0.43           C
    /// ATOM     21  CA  THR A   2      -9.669  -0.447   4.998  1.00  0.19           C
    /// ATOM     35  CA  TYR A   3      -7.173  -2.314   2.811  1.00  0.08           C
    /// ATOM     56  CA  LYS A   4      -3.922  -3.881   4.044  1.00  0.10           C
    /// ATOM     78  CA  LEU A   5      -0.651  -2.752   2.466  1.00  0.11           C
    /// ATOM     97  CA  ILE A   6       2.338  -5.105   2.255  1.00  0.13           C
    /// ATOM      2  CA  MET B   1     -13.296   0.028   3.924  1.00  0.43           C
    /// ATOM     21  CA  THR B   2      -9.669  -0.447   4.998  1.00  0.19           C
    /// ATOM     35  CA  TYR B   3      -7.173  -2.314   2.811  1.00  0.08           C
    /// ATOM     56  CA  LYS B   4      -3.922  -3.881   4.044  1.00  0.10           C
    /// ATOM     78  CA  LEU B   5      -0.651  -2.752   2.466  1.00  0.11           C";
    /// let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_txt.as_bytes())).unwrap();
    /// let strctr = deposit.structure().unwrap();
    /// let first = ResidueId::new("A", 4, ' ');
    /// let last = ResidueId::new("B", 2, ' ');
    /// let mut iterator = strctr.atoms_in_range(first, last);
    /// assert_eq!(iterator.count(), 5);
    /// ```
    pub fn atoms_in_range(&self, first_res: ResidueId, last_res: ResidueId) -> impl Iterator<Item = &PdbAtom> {
        let check = ByResidueRange::new(first_res, last_res);
        self.atoms.iter().filter(move |&a| check.check(a))
    }

    /// Returns ResidueId at a given position in the structure.
    ///
    /// Residue indexes go continuously from 0 through all the chains
    pub fn residue_by_index(&self, index: usize) -> Option<&ResidueId> {
        if index < self.residue_ids.len() {
            Some(&self.residue_ids[index])
        } else {
            None
        }
    }

    /// Iterates over [`ResidueId`](ResidueId)s from a given range of residues
    pub fn residues_in_range<'a>(&'a self, chain_id: &'a str, residue_name: &'a str) -> impl Iterator<Item = &'a ResidueId> + 'a {

        self.residue_ids
            .iter()
            .zip(self.atoms_for_residue_id.iter())
            .filter(move |(res_id, range)| {
                res_id.chain_id == chain_id
                    && range.start < self.atoms.len() // range safety
                    && self.atoms[range.start].res_name == residue_name
            })
            .map(|(res_id, _)| res_id)
    }

    /// Finds all residues in a given chain by their 3-letter name.
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::{Deposit, PDBError, ResidueId};
    /// # fn main() -> Result<(), PDBError> {
    /// let deposit = Deposit::from_file("./tests/test_files/2gb1.cif")?;
    /// let mut strctr = deposit.structure().unwrap();
    /// let n_thr = strctr.residues_in_range("A", "THR").count();
    /// assert_eq!(n_thr, 11);
    /// # Ok(())
    /// # }
    /// ```
    pub fn find_residues_by_name<'a>(&'a self, chain_id: &'a str, residue_name: &'a str) -> impl Iterator<Item=&'a ResidueId> + 'a {
        self.residues()
            .iter()
            .filter(move |res_id| res_id.chain_id == chain_id)
            .filter(move |res_id| self.residue_type(res_id).expect("").code3 == residue_name)
    }

    /// Immutable access to a vector of [`ResidueId`](ResidueId) object for each residue of this [`Structure`](Structure)
    ///
    /// # Examples
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// use bioshell_seq::chemical::StandardResidueType;
    /// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    ///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let strctr = Structure::from_atoms("1xyz", atoms);
    /// assert_eq!(strctr.residues().len(), 2);
    /// assert_eq!(strctr.residues()[0].res_seq, 68);
    /// assert_eq!(strctr.residues()[1].res_seq, 69);
    /// ```
    pub fn residues(&self) -> &Vec<ResidueId> { &self.residue_ids }

    /// Returns a vector of indexes pointing to [`ResidueId`](ResidueId)s  of a given chain.
    ///
    /// Returned indexes simply comprise the longest entity of that chain; it's
    /// assumed that is the polymer chain. In CIF files water molecules and ligands
    /// are placed into separate entities. In the PDB format they are separated from the polymer entity
    /// by a `TER` record and are marked appropriately by bioshell while loading.
    fn residues_in_polymer(&self, chain_id: &str) -> Vec<usize> {

        let mut counts: HashMap<String, Vec<usize>> = HashMap::new();
        for (i_res, res_id) in self.residue_ids.iter().enumerate() {
            if res_id.chain_id != chain_id { continue }
            let first_atom: &PdbAtom = &self.atoms[self.atoms_for_residue_id[i_res].start];
            if first_atom.res_name == "HOH" { continue }
            counts.entry(first_atom.entity_id.clone()).or_insert(vec![]).push(i_res);
        }
        let mut max_len = 0;
        let mut max_id: Option<String> = None;
        for (id, v) in counts.iter() {
            if v.len() > max_len {
                max_len = v.len();
                max_id = Some(id.clone());
            }
        }
        return counts.remove(max_id.as_ref().unwrap()).unwrap();
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue of a given chain
    ///
    pub fn residue_in_chain(&self, chain_id: &str) -> Vec<ResidueId> {

        Structure::residue_ids_from_atoms(self.atoms.iter().filter(|&a| a.chain_id==chain_id))
    }

    /// Returns the chemical type of residue as a [`ResidueType`] object.
    ///
    /// Results in an [`PDBError`] if the type of the residue hasn't been registered in the [`ResidueTypeManager`].
    ///
    /// # Example
    /// ```
    /// # use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// use bioshell_seq::chemical::StandardResidueType;
    /// let pdb_lines = vec!["ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C"];
    /// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let strctr = Structure::from_atoms("1xyz", atoms);
    /// let res_type = strctr.residue_type(&ResidueId::new("A", 68, ' ')).unwrap();
    /// assert_eq!(res_type.code3, "ALA");
    /// assert_eq!(res_type.parent_type, StandardResidueType::ALA);
    /// ```
    pub fn residue_type(&self, res_id: &ResidueId) -> Result<ResidueType, PDBError> {

        // --- check if such a residue has at least one atom in this struct
        if let Some(atom) = self.atoms().iter().find(|&a| res_id.check(a)) {
            if let Some(res_type) = KNOWN_RESIDUE_TYPES.lock().unwrap().by_code3(&atom.res_name) {
                return Ok(res_type.clone());    // --- atom exists and its residue type is known
            } else {                            // --- atom exists but its residue type is NOT known
                return Err(PDBError::UnknownResidueType { res_type: atom.res_name.clone()});
            }
        } else {                        // --- atom doesn't exist
            return Err(NoSuchResidue { res_id: res_id.clone() });
        }
    }

    /// Returns the secondary structure a given residue.
    ///
    /// This method returns [`SecondaryStructureTypes`](SecondaryStructureTypes) enum variant
    /// to define the type of secondary structure element a given residue belongs to. The enum stores
    /// also the index of a secondary structure element in the structure.
    /// Use [`Structure::secondary()`] method to get the secondary structure of a full chain.
    ///
    /// # Examples
    /// Check the secondary structure of a residue in a PDB deposit:
    /// ```
    /// # use bioshell_pdb::{Deposit, PDBError, ResidueId, SecondaryStructureTypes};
    /// # fn main() -> Result<(), PDBError> {
    /// # let deposit = Deposit::from_file("./tests/test_files/2gb1.cif")?;
    /// let mut strctr = deposit.structure().unwrap();
    /// let maybe_helix = strctr.residue_secondary(&ResidueId::new("A", 27, ' '))?;
    /// assert!(matches!(maybe_helix, SecondaryStructureTypes::RightAlphaHelix(_)));
    /// let maybe_strand = strctr.residue_secondary(&ResidueId::try_from("A:4")?)?;
    /// assert!(matches!(maybe_strand, SecondaryStructureTypes::Strand(_)));
    /// # Ok(())
    /// # }
    /// ```
    ///
    /// The example below counts how many times `'H'`, `'E'` and `'C'` appears in the secondary
    /// structure of the 2gb1 deposit, loaded from a local file:
    /// ```
    /// # use bioshell_pdb::{Deposit, PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// use std::collections::HashMap;
    /// let deposit = Deposit::from_file("./tests/test_files/2gb1.cif")?;
    /// let mut strctr = deposit.structure().unwrap();
    /// let counts = strctr.residues().iter()
    ///     .map(|rid| strctr.residue_secondary(rid).expect("ERROR").hec_code())
    ///     .fold(HashMap::new(), |mut acc, item| {
    ///         *acc.entry(item).or_insert(0) += 1;
    ///         acc
    ///     });
    /// assert_eq!(counts[&b'H'], 15);
    /// assert_eq!(counts[&b'E'], 22);
    /// assert_eq!(counts[&b'C'], 19);
    /// # Ok(())
    /// }
    /// ```
    pub fn residue_secondary(&self, res_id: &ResidueId) -> Result<SecondaryStructureTypes, PDBError> {

        // --- check if such a residue has at least one atom in this struct
        if let Some(atom) = self.atoms().iter().find(|&a| res_id.check(a)) {
            return Ok(atom.secondary_struct_type.clone());
        } else {                        // --- atom doesn't exist
            return Err(NoSuchResidue { res_id: res_id.clone() });
        }
    }
    /// Provides a sequence of a given chain.
    ///
    /// The sequence contains only in the residues found in atoms of this structure; some of its residues
    /// may be missing. The original sequence of a chain can be obtained from the respective [`Entity`](crate::Entity) object.
    ///
    /// This method includes only the residues included in a [`Polymer`](crate::EntityType::Polymer)
    /// entity of the requested chain. Because such an information is not provided by the PDB file format,
    /// in that case the method includes residues listed in the requested chain until the `TER` record.
    /// ```
    /// # use std::io::BufReader;
    /// # use bioshell_pdb::Deposit;
    /// # use bioshell_pdb::PDBError;
    /// # fn main() -> Result<(), PDBError> {
    /// let pdb_data = "
    /// ATOM   4603  CA  GLN H 244      32.033  12.635  17.439  1.00 50.48           C
    /// ATOM   4620  CA  PHE H 245      30.578  11.461  14.097  1.00 53.40           C
    /// ATOM   4640  CA  GLY H 246      27.584  13.778  13.222  0.81 63.93           C
    /// ATOM   4647  CA  GLU H 247      23.988  14.633  14.349  0.41 70.87           C
    /// TER    4662      GLU H 247
    /// HETATM 4835 NA    NA H 409       5.622 -14.085  31.049  1.00 32.84          NA
    /// HETATM 4836 CA    CA H 521      18.674 -16.434  38.353  0.36  7.57          CA
    /// ";
    /// let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_data.as_bytes()))?;
    /// let strctr = deposit.structure().unwrap();
    /// let seq = strctr.sequence("H");
    /// assert_eq!("QFGE", &seq.to_string(10));
    /// # Ok(())
    /// # }
    /// ```
    pub fn sequence(&self, chain_id: &str) -> Sequence {

        let res_ids = self.residues_in_polymer(chain_id);
        let rtm = ResidueTypeManager::get();
        let mut residue_sequence: Vec<u8> = vec![];
        for i_res in res_ids {
            let first_atom: &PdbAtom = &self.atoms[self.atoms_for_residue_id[i_res].start];
            // --- if the monomer type of this residue has been already registered,
            // --- use its code1, otherwise use 'X'.
            let mut code_1 = b'X';
            if let Some(res_type) = rtm.by_code3(&first_atom.res_name) {
                code_1 = res_type.parent_type.code1() as u8;
            }
            residue_sequence.push(code_1);
        }

        return Sequence::from_attrs(format!("{}:{}", &self.id_code, chain_id), residue_sequence);
    }

    /// Provides a secondary structure of a given chain.
    ///
    /// This method includes only the residues included in a [`Polymer`](Polymer) entity of the requested chain.
    /// Because such an information is not provided by the PDB file format, in that case the method
    /// includes residues listed in the requested chain until the `TER` record.
    ///
    /// You can use [`Structure::residue_secondary()`] method to get the secondary structure of a single residue.
    pub fn secondary(&self, chain_id: &str) -> SecondaryStructure {
        let mut sec: Vec<SecondaryStructureTypes> = vec![];
        let res_ids = self.residues_in_polymer(chain_id);

        for i_res in res_ids {
            let first_atom: &PdbAtom = &self.atoms[self.atoms_for_residue_id[i_res].start];
            sec.push(first_atom.secondary_struct_type.clone());
        }

        return SecondaryStructure::from_attrs(sec);
    }

    /// Removes ligands and waters from a structure.
    ///
    /// This method assumes a polymer molecule is listed first in a PDB chain, ligands and water molecules
    /// come after it. This method removes all the entities other than the polymer chain entity.
    /// ```
    /// # use bioshell_pdb::{Deposit, PdbAtom, ResidueId, Structure};
    /// # use std::io::BufReader;
    /// let pdb_lines ="
    /// ATOM   4378  CA  HIS D 146      14.229  -1.501  26.683  1.00 31.89           C
    /// TER    4388      HIS D 146
    /// HETATM 4562 FE   HEM D 148      -1.727   4.699  23.942  1.00 15.46          FE;
    /// ATOM   4378  CA  HIS E 146      14.229  -1.501  26.683  1.00 31.89           C
    /// TER    4388      HIS E 146
    /// HETATM 4562 FE   HEM E 148      -1.727   4.699  23.942  1.00 15.46          FE";
    ///
    /// let deposit = Deposit::from_pdb_reader(BufReader::new(pdb_lines.as_bytes())).unwrap();
    /// let mut strctr = deposit.structure().unwrap();
    /// strctr.remove_ligands();
    /// assert_eq!(strctr.count_atoms(), 2);
    /// ```
    pub fn remove_ligands(&mut self) {
        let mut first_id_map: HashMap<String, String> = HashMap::new();

        self.atoms
            .retain(|atom| {
                // Check if we've seen this entity before
                if let Some(first_id) = first_id_map.get(&atom.chain_id) {
                    // If we have, only keep it if the id matches the first one seen
                    &atom.entity_id == first_id
                } else {
                    // If we haven't seen this entity, store the id and keep the element
                    first_id_map.insert(atom.chain_id.clone(), atom.entity_id.clone());
                    true
                }
            });
    }

    /// Creates a vector of [`ResidueId`](ResidueId) object for each residue found in a given vector of atoms
    ///
    /// ```
    /// # use bioshell_pdb::{PdbAtom, Structure};
    /// let pdb_lines = vec!["ATOM    514  N   ALA A  68      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  68      25.790  28.757  29.513  1.00 16.12           C",
    ///                      "ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N",
    ///                      "ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"];
    /// let atoms: Vec<PdbAtom> = pdb_lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
    /// let res_idx = Structure::residue_ids_from_atoms(atoms.iter());
    /// assert_eq!(format!("{}", res_idx[0]), "A:68 ");
    /// assert_eq!(format!("{}", res_idx[1]), "A:69 ");
    /// assert_eq!(res_idx.len(), 2);
    /// ```
    pub fn residue_ids_from_atoms<'a>(atoms: impl Iterator<Item = &'a PdbAtom>) -> Vec<ResidueId> {
        // --- predicate used to check whether we are entering a new residue
        let same_res = SameResidue{};
        // --- turn a given iterator into a peekable one
        let mut peek_iter = atoms.peekable();
        // --- take the first PdbAtom, if there is any
        let maybe_first = peek_iter.peek();

        if let Some(&first) = maybe_first {
            // --- take the ResidueIdx of the first PdbAtom
            let first_idx = ResidueId::try_from(first).unwrap();
            let mut ret: Vec<ResidueId> = peek_iter.tuple_windows()
                .filter(|(a,b)| !same_res.check(&a, &b))
                .map(|(_, b)| ResidueId::try_from(b).unwrap()).collect();

            // --- insert the first ResidueIdx which we actually juped over during the iteration above
            ret.insert(0, first_idx);

            return ret;
        }
        // --- return empty vector if there are no atoms in the given iterator
        return Vec::new();
    }

    /// Sorts all the atoms of this structure.
    pub fn sort(&mut self) { self.atoms.sort(); }

    /// Returns true if all the atoms of this structure are sorted
    pub fn is_sorted(&self) -> bool {
        self.atoms.windows(2).all(|w| w[0] <= w[1])
    }

    /// Returns the index of a residue given its [`ResidueId`](ResidueId)
    ///
    pub(crate) fn residue_pos(&self, which_res: &ResidueId) -> Result<usize, PDBError> {
        // --- check if residue_ids are sorted
        if let Ok(pos) = self.residue_ids.binary_search(which_res) {
            return Ok(pos);
        }
        // --- otherwise go through the vector of residue ids and find the position
        match self.residue_ids.iter().position(|r| r == which_res) {
            Some(pos) => Ok(pos),
            None => Err(NoSuchResidue { res_id: which_res.clone() })
        }
    }

    /// Updates the internal structure of the struct
    ///
    /// This method should be called after any change to the atoms of this Structure
    pub(crate) fn update(&mut self) {
        if self.atoms.is_empty() { return; }

        // --- sort atoms, just in case
        // self.sort();
        // --- assign atom range for every residue
        self.setup_atom_ranges();
    }

    pub(crate) fn setup_atom_ranges(&mut self) {
        self.atoms_for_residue_id.clear();
        self.residue_ids.clear();
        let mut first_of_res: usize = 0;
        let mut res_id = ResidueId::try_from(&self.atoms[0]).unwrap();  // never fails
        for i_atom in 1..self.atoms.len() {
            if ! same_residue_atoms(&self.atoms[first_of_res], &self.atoms[i_atom]) {
                self.atoms_for_residue_id.push(first_of_res..i_atom);
                self.residue_ids.push(res_id);
                first_of_res = i_atom;
                res_id = ResidueId::try_from(&self.atoms[first_of_res]).unwrap();
            }
        }
        self.atoms_for_residue_id.push(first_of_res..self.atoms.len());
        self.residue_ids.push(res_id);
    }

    #[allow(dead_code)]
    /// Borrows the very last atom of this [`Structure`]
    fn last_atom(&self) -> &PdbAtom { &self.atoms[self.atoms.len()-1] }

}

/// Write a given structure in the PDB format.
///
/// The structure stored in a file may differ from the given object, as it may need to be adapted to the PDB file format.
/// For example, if:
///  - a given structure has a chain with its name longer than a single character, all chains will be renamed
pub fn write_pdb(strctr: &Structure, mut outstream: Box<dyn Write>) -> bool {

    let mut if_rename_chains = false;
    for chain_id  in &strctr.chain_ids() {
        if chain_id.len() > 1 {
            if_rename_chains = true;
            break
        }
    }
    let new_chain_codes: Vec<char> = ('A'..='Z').chain('0'..='9').collect();
    for (i, chain_id)  in strctr.chain_ids().iter().enumerate() {
        for atom in strctr.atoms_in_chain(chain_id) {
            let mut a = atom.clone();
            if if_rename_chains {
                a.chain_id = new_chain_codes[i].to_string();
            }
            writeln!(outstream, "{}", a).unwrap();
        }
    }

    outstream.flush().unwrap();

    return true;
}


