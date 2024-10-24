use crate::{ResidueId, Structure};
use crate::monomers::MonomerManager;
use bioshell_seq::chemical::StandardResidueType::{TYR, PHE, TRP, HIS};

/// A handy filter to process residues of a [`Structure`](crate::Structure) with iterators.
///
/// Structs implementing [`ResidueFilter`](crate::residue_filters::ResidueFilter) trait can be used as predicates
/// while filtering Rust iterators. Example below shows how to iterate over aromatic amino acid residues
/// of a protein chain:
/// ```
/// # use bioshell_pdb::{PDBError, Deposit};
/// # fn main() -> Result<(), PDBError> {
/// use bioshell_pdb::residue_filters::{IsAromaticAA, ResidueFilter};
/// # let cif_data = include_str!("../tests/test_files/2fdo.cif");
/// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
/// let strctr = deposit.structure();
/// let n_aro = strctr.residue_ids().iter().filter(|ri| IsAromaticAA.check(&strctr, &ri)).count();
/// assert_eq!(n_aro, 28);
/// # Ok(())
/// # }
/// ```
pub trait ResidueFilter {
    /// Returns `true` if this predicate is satisfied
    fn check(&self, strctr: &Structure, ri: &ResidueId) -> bool;
}

/// Returns `true` if the residue has all backbone atoms.
pub struct IsAromaticAA;

impl ResidueFilter for IsAromaticAA {
    fn check(&self, strctr: &Structure, ri: &ResidueId) -> bool {
        if let Ok(monomer) = &strctr.residue_type(ri) {
            return matches!(monomer.parent_type, TYR | TRP | PHE | HIS);
        }
        return false
    }
}

const BB_ATOMS: [&str; 4] = [" N  ", " CA ", " C  ", " O  "];

/// Returns `true` if the residue has all backbone atoms
pub struct HasCompleteBackbone;

impl ResidueFilter for HasCompleteBackbone {
    /// Creates a new filter that accepts residues with complete backbone atoms
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::residue_filters::{HasCompleteBackbone, ResidueFilter};
    /// use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    /// let filter = HasCompleteBackbone;
    /// # let mut strctr = Structure::new("1xyz");
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    514  N   ALA A  69      26.532  28.200  28.365  1.00 17.85           N"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    515  CA  ALA A  69      25.790  28.757  29.513  1.00 16.12           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    516  C   ALA A  69      26.891  29.054  30.649  1.00 15.28           C"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    517  O   ALA A  69      26.657  29.867  31.341  1.00 20.90           O"));
    /// # strctr.push_atom(PdbAtom::from_atom_line("ATOM    518  CB  ALA A  69      25.155  27.554  29.987  1.00 21.91           C"));
    /// assert!(filter.check(&strctr, &ResidueId::new("A", 69, ' ')));
    /// ```
    fn check(&self, strctr: &Structure, ri: &ResidueId) -> bool {
        for atom_name in BB_ATOMS {
            if strctr.atom(ri, atom_name).is_err() { return false}
        }
        return true;
    }
}

/// Returns `true` if the residue has all its atoms, **including hydrogens**.
///
/// Leaving atoms (those removed when a polymer link was formed) are also not counted.
pub struct HasAllAtoms;
impl ResidueFilter for HasAllAtoms {

    fn check(&self, strctr: &Structure, ri: &ResidueId) -> bool {
        let mgr = MonomerManager::get();
        let atoms = strctr.atoms_in_residue(ri).unwrap();
        return if let Some(first_atom) = atoms.peekable().peek() {
            let residue_def = mgr.by_code3(&first_atom.res_name);
            if residue_def.is_none() { return false; }
            let residue_def = residue_def.unwrap();
            let atoms_found = strctr.atoms_in_residue(ri).iter().count();
            atoms_found >= residue_def.count_residue_atoms()
        } else { false }
    }
}

/// Returns `true` if the residue has all its heavy atoms, i.e. hydrogens are not counted.
///
/// Returns `false` when :
///     - the residue referred by the given ``ResidueId`` doesn't exist in a structure
///     - the residue type is not registered in the monomer manager
///     - the residue has too few atoms
pub struct HasAllHeavyAtoms;
impl ResidueFilter for HasAllHeavyAtoms {

    fn check(&self, strctr: &Structure, ri: &ResidueId) -> bool {

        // --- check if we have the residue and its type
        let res_type = strctr.residue_type(ri);
        if res_type.is_err() { return false; }
        let res_type = res_type.unwrap();

        // --- check if we have structure definition for that residue type
        let mgr = MonomerManager::get();

        let res_def = mgr.by_code3(&res_type.code3);
        if res_def.is_none() { return false; }
        let res_def = res_def.unwrap();

        // --- list atoms of the tested residue
        let atoms = strctr.atoms_in_residue(ri);
        if atoms.is_err() { return false; }
        let atoms = atoms.unwrap();

        let mut atom_cnt = 0;
        let h: Option<String> = Some("H".to_string());
        for atom in atoms {
            if atom.element.is_none() { continue }
            if atom.element != h { atom_cnt += 1; }
        }

        return atom_cnt >= res_def.count_residue_heavy();
    }
}

/// Tests a condition that relates two residues of a [`Structure`](crate::Structure)
///
pub trait ResidueFilter2 {
    /// Returns `true` if this predicate is satisfied
    fn check(&self, strctr: &Structure, ri: &ResidueId, rj: &ResidueId) -> bool;
}

/// Returns `true` if the two residues are connected by a peptide bond.
///
/// The peptide bond if detected when the ``first`` residue has the carbonyl ``C``,
/// the ``second`` residue has its amide ``N`` atom, and the distance between them is less than 1.4 Å.
pub struct ArePeptideBonded;

impl ResidueFilter2 for ArePeptideBonded {
    fn check(&self, strctr: &Structure, first: &ResidueId, second: &ResidueId) -> bool {
        if let Ok(c)  = strctr.atom(first, " C  ") {
            if let Ok(n)  = strctr.atom(second, " N  ") {
                if c.pos.distance_to(&n.pos) < 1.4 { return true; }
            }
        }  return false;
    }
}
