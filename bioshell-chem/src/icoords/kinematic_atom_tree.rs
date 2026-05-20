use std::collections::{HashSet, VecDeque};
use std::ops::Index;

use bioshell_core::{dihedral_angle4, planar_angle3, Vec3};
use crate::{ChemErrors, Molecule};
use crate::ChemErrors::IncorrectNumberOfAtoms;

/// Numerical internal coordinates for a single atom.
///
/// The fields correspond to distance `d`, planar angle `alpha`, and dihedral
/// angle `phi`. Values are determined by the caller; distance follows
/// the Cartesian coordinate units in Angstroms, while angles are stored in radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InternalCoordinate {
    pub d: f64,
    pub alpha: f64,
    pub phi: f64,
}

impl Default for InternalCoordinate {
    fn default() -> Self {
        Self {
            d: 1.65,
            alpha: 108.0_f64.to_radians(),
            phi: 180.0_f64.to_radians(),
        }
    }
}

/// A reference frame to locate an atom by three other reference atoms.
///
/// The three atom, indexed by `a`, `b` and `c` will be used as a reference frame for a given `atom`.
/// The position of the `atom` will be defined by three internal coordinates:
///    * `r` - distance between `c` and the `atom`
///    * `planar` - planar angle between `b`, `c` and `atom`
///    * `dihedral` - dihedral angle between `a`, `b`, `c` and `atom`
/// ```text
///         c----atom
///        /
///       /
/// a----b
/// ```
///
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KinematicAtom {
    pub atom: usize,
    pub a: usize,
    pub b: usize,
    pub c: usize,
}


/// A molecular internal-coordinate tree.
///
/// The chain stores the topological definition of each atom of a molecule
/// in relation to the other atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct KinematicAtomTree {
    atoms: Vec<KinematicAtom>,
}

impl KinematicAtomTree {
    /// Returns the number of atoms in the chain.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Returns `true` if the chain contains no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Returns a read-only slice of all kinematic atoms.
    pub fn atoms(&self) -> &[KinematicAtom] {
        &self.atoms
    }

    /// Returns a read-only reference to the kinematic atom at `index`.
    pub fn atom(&self, index: usize) -> &KinematicAtom {
        &self.atoms[index]
    }


    /// Creates a kinematic atom chain from a molecular graph.
    ///
    /// The first three atoms indexed by `i`, `j`, and `k` are provided explicitly by the caller.
    /// This allows chemically meaningful initialization of the internal-coordinate chain.
    /// Remaining atoms are added by graph traversal, with heavy atoms preferred
    /// over hydrogens when several traversal/reference choices are possible.
    ///
    /// ```
    /// # use bioshell_chem::{ChemErrors, Molecule, BondType, Atom, Element};
    /// # use bioshell_chem::icoords::KinematicAtomTree;
    /// # fn main() -> Result<(), ChemErrors> {
    /// # let mut benzene = Molecule::new("benzene");
    /// # for i in 0..6 { benzene.add_atom(Atom::neutral(i, Element::C))?; }
    /// # for i in 0..6 { benzene.bind_atoms(i, (i + 1) % 6, BondType::Aromatic)?; }
    ///
    /// let chain = KinematicAtomTree::from_molecule(&benzene, 0, 1, 2)?;
    /// assert_eq!(chain.len(), 6);
    /// assert_eq!(chain.atom(0).atom, 0);
    /// assert_eq!(chain.atom(1).atom, 1);
    /// assert_eq!(chain.atom(2).atom, 2);
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_molecule(mol: &Molecule, i: usize, j: usize, k: usize) -> Result<KinematicAtomTree, ChemErrors> {
        mol.get_atom(i).ok_or(ChemErrors::InvalidAtomIndex(i))?;
        mol.get_atom(j).ok_or(ChemErrors::InvalidAtomIndex(j))?;
        mol.get_atom(k).ok_or(ChemErrors::InvalidAtomIndex(k))?;

        if i == j || i == k || j == k {
            return Err(ChemErrors::InvalidReferenceAtoms(i, j, k));
        }

        let mut atoms = Vec::new();
        let mut placed = HashSet::new();
        let mut queue = VecDeque::new();

        atoms.push(KinematicAtom { atom: i, a: i, b: i, c: i });
        atoms.push(KinematicAtom { atom: j, a: i, b: i, c: i });
        atoms.push(KinematicAtom { atom: k, a: i, b: i, c: j });

        placed.insert(i);
        placed.insert(j);
        placed.insert(k);

        queue.push_back(i);
        queue.push_back(j);
        queue.push_back(k);

        while let Some(center) = queue.pop_back() {
            let neighbors = KinematicAtomTree::prioritized_atoms(mol, mol.neighbor_indices(center))?;

            for atom in neighbors {
                if placed.contains(&atom) { continue; }

                let bond_ref = center;
                let angle_ref = KinematicAtomTree::choose_angle_reference(mol, atom, bond_ref, &placed)?;
                let dihedral_ref = KinematicAtomTree::choose_dihedral_reference(mol, atom, bond_ref, angle_ref, &placed)?;

                atoms.push(KinematicAtom { atom, a: dihedral_ref, b: angle_ref, c: bond_ref});
                placed.insert(atom);
                queue.push_back(atom);
            }
        }

        if atoms.len() != mol.count_atoms() {
            return Err(ChemErrors::InvalidAtomIndex(atoms.len()));
        }

        Ok(KinematicAtomTree { atoms })
    }

    pub fn get_icoords(&self, pos: &[Vec3]) -> Result<Vec<InternalCoordinate>, ChemErrors> {
        let mut out = vec![InternalCoordinate::default(); self.len()];
        self.get_icoords_into(pos, &mut out)?;

        return Ok(out);
    }

    pub fn get_icoords_into(&self, pos: &[Vec3], out: &mut [InternalCoordinate]) -> Result<(), ChemErrors> {
        if out.len() != self.len() {
            return Err(IncorrectNumberOfAtoms(self.len(), out.len()));
        }
        if pos.len() != self.len() {
            return Err(IncorrectNumberOfAtoms(self.len(), pos.len()));
        }

        out[0] = InternalCoordinate { d: 0.0, alpha: 0.0, phi: 0.0 };
        let ka = &self.atoms[1];
        out[1] = InternalCoordinate { d: pos[ka.atom].distance_to(&pos[ka.c]), alpha: 0.0, phi: 0.0};
        let ka = &self.atoms[2];
        out[2] = InternalCoordinate {
            d: pos[ka.atom].distance_to(&pos[ka.c]),
            alpha: planar_angle3(&pos[ka.b], &pos[ka.c], &pos[ka.atom]),
            phi: 0.0,
        };
        for row in 3..self.len() {
            let ka = &self.atoms[row];
            out[row] = InternalCoordinate {
                d: pos[ka.atom].distance_to(&pos[ka.c]),
                alpha: planar_angle3(&pos[ka.b], &pos[ka.c], &pos[ka.atom]),
                phi: dihedral_angle4(&pos[ka.a], &pos[ka.b], &pos[ka.c], &pos[ka.atom]),
            };
        }
        return Ok(());
    }

    fn atom_priority(mol: &Molecule, atom_idx: usize) -> Result<(usize, usize, usize), ChemErrors> {

        let center = mol.get_atom(atom_idx).ok_or_else(|| ChemErrors::InvalidAtomIndex(atom_idx))?;
        let heavy_rank = if center.if_hydrogen() { 1 } else { 0 };
        let degree_rank = usize::MAX - mol.neighbor_indices(atom_idx).count();

        Ok((heavy_rank, degree_rank, atom_idx))
    }

    fn prioritized_atoms(mol: &Molecule, atoms: impl Iterator<Item = usize>) -> Result<Vec<usize>, ChemErrors>  {
        let mut out: Vec<(usize, (usize, usize, usize))> = atoms
            .map(|idx| Ok((idx, KinematicAtomTree::atom_priority(mol, idx)?)))
            .collect::<Result<_, ChemErrors>>()?;

        out.sort_by_key(|&(_, priority)| priority);

        return Ok(out.into_iter().map(|(idx, _)| idx).collect());
    }

    fn choose_angle_reference(
        mol: &Molecule,
        atom: usize,
        bond_ref: usize,
        placed: &HashSet<usize>,
    ) -> Result<usize, ChemErrors> {
        let candidates = KinematicAtomTree::prioritized_atoms(
            mol,
            mol.neighbor_indices(bond_ref)
                .filter(|&idx| idx != atom && placed.contains(&idx)),
        )?;

        Ok(candidates.first().map(|&idx| idx).unwrap_or(bond_ref))
    }

    fn choose_dihedral_reference(
        mol: &Molecule,
        atom: usize,
        bond_ref: usize,
        angle_ref: usize,
        placed: &HashSet<usize>,
    ) -> Result<usize, ChemErrors> {
        let candidates = KinematicAtomTree::prioritized_atoms(
            mol,
            mol.neighbor_indices(angle_ref)
                .filter(|&idx| idx != atom && idx != bond_ref && placed.contains(&idx)),
        )?;

        Ok(candidates.first().copied().unwrap_or(angle_ref))
    }
}


impl Index<usize> for KinematicAtomTree {
    type Output = KinematicAtom;

    fn index(&self, index: usize) -> &Self::Output {
        &self.atoms[index]
    }
}