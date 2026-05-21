use std::fmt;
use std::fmt::Display;

use bioshell_core::Vec3;
use crate::icoords::{InternalCoordinate, KinematicAtom, KinematicAtomTree};
use crate::{ChemErrors, Molecule};

pub struct ZMatrix<'a> {
    molecule: &'a mut Molecule,
    kinematic_tree: KinematicAtomTree,
    icoords: Vec<InternalCoordinate>,
}

impl<'a> ZMatrix<'a> {
    pub fn from_molecule(molecule: &'a mut Molecule, i: usize, j: usize, k: usize) -> Result<Self, ChemErrors> {
        let chain = KinematicAtomTree::from_molecule(molecule, i, j, k)?;

        let cartesian: &Vec<Vec3> = molecule.positions();

        let icoords = chain.get_icoords(cartesian)?;

        return Ok(Self { molecule, kinematic_tree: chain, icoords});
    }

    /// Number of atoms  in the Z-matrix.
    ///
    /// The returned value is the same as the number of atoms in the underlying Molecule.
    pub fn len(&self) -> usize {
        self.kinematic_tree.len()
    }

    pub fn atom_definition(&self, atom_index: usize) -> &KinematicAtom {
        &self.kinematic_tree[atom_index]
    }

    // Internal coordinates of one atom
    pub fn internal_coordinate(&self, atom_index: usize) -> &InternalCoordinate {
        &self.icoords[atom_index]
    }

    // Cartesian coordinates of one atom, read from the underlying Molecule
    pub fn cartesian_coordinate(&self, atom_idx: usize) -> &Vec3 {
        self.molecule.pos(atom_idx)
    }

    pub fn internal_coordinates(&self) -> &[InternalCoordinate] {
        &self.icoords
    }

    pub fn atom_definitions(&self) -> &KinematicAtomTree {
        &self.kinematic_tree
    }

    /// Cartesian coordinates of all the atoms from the underlying Molecule
    pub fn cartesian_coordinates(&self, ) -> &Vec<Vec3> {
        self.molecule.positions()
    }

    pub fn write<W: fmt::Write>(&self, mut writer: W) -> Result<(), fmt::Error> {

        for i in 0..self.icoords.len() {
            let ipos = &self.icoords[i];
            let ka = &self.kinematic_tree[i];
            let atom = self.molecule.get_atom(ka.atom).ok().expect("Inconsistent Z-matrix: atom index out of bounds");

            writeln!(
                writer,
                "{:3} {:2} {:3} {:3} {:3} : {:5.3} {:6.2} {:7.2}",
                ka.atom,
                atom.element(),
                ka.a,
                ka.b,
                ka.c,
                ipos.d,
                ipos.alpha.to_degrees(),
                ipos.phi.to_degrees()
            )?;
        }

        Ok(())
    }
}

impl Display for ZMatrix<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write(f).map_err(|_| fmt::Error)
    }
}