use std::fmt;
use std::fmt::Display;

use bioshell_core::Vec3;
use crate::icoords::{InternalCoordinate, KinematicAtom, KinematicAtomTree};
use crate::{ChemErrors, Molecule};
use crate::icoords::nerf::{create_stub, restore_atom};

pub struct ZMatrix<'a> {
    molecule: &'a mut Molecule,
    kinematic_tree: KinematicAtomTree,
    icoords: Vec<InternalCoordinate>,
}

impl<'a> ZMatrix<'a> {
    pub fn from_molecule(molecule: &'a mut Molecule, i: usize, j: usize, k: usize) -> Result<Self, ChemErrors> {

        let kat = KinematicAtomTree::from_molecule(molecule, i, j, k)?;
        let cartesian: &Vec<Vec3> = molecule.positions();
        let icoords = kat.get_icoords(cartesian)?;

        return Ok(Self { molecule, kinematic_tree: kat, icoords});
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

    /// Updates the internal coordinates of this z-matrix and recalculates the Cartesian coordinates.
    ///
    /// # Example
    /// ```
    /// # use bioshell_chem::ChemErrors;
    /// # fn main() -> Result<(), ChemErrors> {
    /// # use bioshell_chem::icoords::ZMatrix;
    /// # use bioshell_chem::load_molecule;
    /// let mut EtOH = load_molecule("./tests/test_files/EOH.cif")?;
    /// let mut zmatrix = ZMatrix::from_molecule(&mut EtOH, 0, 1, 2)?;
    /// let new_ic = zmatrix.internal_coordinates().to_vec();
    /// zmatrix.set_internals(&new_ic)?;
    /// for a in EtOH.positions() { eprintln!("{}", &a) }
    ///
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_internals(&mut self, internal_coordinates: &[InternalCoordinate]) -> Result<(), ChemErrors> {

        for i in 0..internal_coordinates.len() {
            self.icoords[i].d = internal_coordinates[i].d;
            self.icoords[i].alpha = internal_coordinates[i].alpha;
            self.icoords[i].phi = internal_coordinates[i].phi;
        }

        self.recalculate_internals()?;

        Ok(())
    }

    pub fn atom_definitions(&self) -> &KinematicAtomTree {
        &self.kinematic_tree
    }

    /// Cartesian coordinates of all the atoms from the underlying Molecule
    pub fn cartesian_coordinates(&self, ) -> &Vec<Vec3> {
        self.molecule.positions()
    }

    /// Updates the Cartesian coordinates of all atoms in the underlying Molecule and recalculates the internal coordinates.
    ///
    /// # Example
    /// ```
    /// # use bioshell_chem::ChemErrors;
    /// # fn main() -> Result<(), ChemErrors> {
    /// # use bioshell_chem::icoords::ZMatrix;
    /// # use bioshell_chem::load_molecule;
    /// let mut EtOH1 = load_molecule("./tests/test_files/EOH.cif")?;
    /// let mut zmatrix = ZMatrix::from_molecule(&mut EtOH1, 0, 1, 2)?;
    /// # let new_cartesian = zmatrix.cartesian_coordinates().clone();
    /// zmatrix.set_cartesians(&new_cartesian)?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn set_cartesians(&mut self, cartesian: &[Vec3]) -> Result<(), ChemErrors> {
        if cartesian.len() != self.molecule.count_atoms() {
            return Err(ChemErrors::IncorrectNumberOfAtoms(cartesian.len(), self.molecule.count_atoms()));
        }

        for i in 0..cartesian.len() {
            self.molecule.set_pos3(i, cartesian[i].x, cartesian[i].y, cartesian[i].z);
        }

        self.kinematic_tree.get_icoords_into(cartesian, &mut self.icoords)?;

        Ok(())
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

    fn recalculate_internals(&mut self) -> Result<(), ChemErrors> {

        // --- make an alias to shorten the code
        let def = &self.kinematic_tree;

        // --- two temporary vectors to store reconstructed positions
        let mut r1 = Vec3::default();
        let mut r2 = Vec3::default();

        // --- create the stub = reconstruct the atoms (2) and (3)
        let a =  &self.molecule.positions()[def[0].atom];
        create_stub(&a, self.icoords[1].d, self.icoords[2].d, self.icoords[2].alpha,
                &mut r1, &mut r2);
        self.molecule.set_pos3(def[1].atom, r1.x, r1.y, r1.z);
        self.molecule.set_pos3(def[2].atom, r2.x, r2.y, r2.z);

        // --- reconstruct all the remaining atoms relative to the stub
        for i in 3..self.icoords.len() {
            let def_i = &def[i];
            let inc_i = &self.icoords[i];
            let a =  &self.molecule.positions()[def_i.a];
            let b =  &self.molecule.positions()[def_i.b];
            let c =  &self.molecule.positions()[def_i.c];

            restore_atom(
                a, b, c,   // Cartesian coordinates of the reference atoms
                inc_i.d, inc_i.alpha, inc_i.phi, // Internal coordinates of the i-th atom
                &mut r1 // Cartesian coordinates of the i-th atom to be updated
            );
            self.molecule.set_pos3(def_i.atom, r1.x, r1.y, r1.z)
        }
        Ok(())
    }
}

impl Display for ZMatrix<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write(f).map_err(|_| fmt::Error)
    }
}