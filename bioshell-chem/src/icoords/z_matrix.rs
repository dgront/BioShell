use bioshell_core::Vec3;
use crate::icoords::{InternalCoordinate, KinematicAtom, KinematicAtomTree};
use crate::{ChemErrors, Molecule};
use crate::ChemErrors::InvalidAtomIndex;

pub struct ZMatrix<'a> {
    molecule: &'a mut Molecule,
    kinematic_tree: KinematicAtomTree,
    icoords: Vec<InternalCoordinate>,
}

impl<'a> ZMatrix<'a> {
    pub fn from_molecule(molecule: &'a mut Molecule, i: usize, j: usize, k: usize) -> Result<Self, ChemErrors> {
        let chain = KinematicAtomTree::from_molecule(molecule, i, j, k)?;

        let cartesian: Vec<Vec3> = molecule
            .atoms()
            .map(|atom| atom.pos().clone())
            .collect();

        let icoords = chain.get_icoords(&cartesian)?;

        return Ok(Self { molecule, kinematic_tree: chain, icoords});
    }

    // Number of atoms / rows in the Z-matrix
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
    pub fn cartesian_coordinate(&self, atom_index: usize) -> Option<Vec3> {
        self.molecule.get_atom(atom_index).map(|atom| atom.pos().clone())
    }

    pub fn internal_coordinates(&self) -> &[InternalCoordinate] {
        &self.icoords
    }

    pub fn atom_definitions(&self) -> &KinematicAtomTree {
        &self.kinematic_tree
    }

    pub fn write<W: std::io::Write>(&self, mut writer: W) -> Result<(), ChemErrors> {

        for i in 0..self.icoords.len() {
            let ipos = &self.icoords[i];
            let ka = &self.kinematic_tree[i];
            let atom = &self.molecule.get_atom(ka.atom)
                .ok_or(InvalidAtomIndex(ka.atom))?;

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