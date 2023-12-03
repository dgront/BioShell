use std::collections::HashMap;
use crate::calc::Vec3;
use crate::nerf::{InternalAtomDefinition, restore_branched_chain};


/// Computes Cartesian coordinates of atoms from theirs internal definition
pub struct KinematicAtomTree {
    r: Vec<f64>,
    planar: Vec<f64>,
    dihedral: Vec<f64>,
    names: HashMap<usize, String>,
    reference_atoms: Vec<[usize;4]>,
    last_atom_by_name: HashMap<String, usize>
}

impl KinematicAtomTree {

    /// Creates new [`KinematicAtomTree`](KinematicAtomTree)  builder.
    pub fn new() -> KinematicAtomTree {
        KinematicAtomTree{ r: vec![], planar: vec![], dihedral: vec![], names: HashMap::default(),
            reference_atoms: vec![], last_atom_by_name: Default::default()
        }
    }

    /// The number of atoms in this tree.
    ///
    /// Returns the number of atoms stored in this kinematic tree
    pub fn len(&self) -> usize { self.reference_atoms.len() }

    /// Name of a given atom.
    ///
    /// Returns the name of an atom as defined in a respective [`InternalAtomDefinition`](InternalAtomDefinition)
    /// object inserted with [`add_atom()`](add_atom()) call
    pub fn name(&self, i: usize) -> Option<&String> { self.names.get(&i) }

    /// Provides the distance used to define `i`-th atom
    pub fn r(&self, i: usize) -> f64 { self.dihedral[i] }

    /// Provides the planar angle used to define `i`-th atom
    pub fn planar(&self, i: usize) -> f64 { self.planar[i] }

    /// Provides the dihedral angle used to define `i`-th atom
    pub fn dihedral(&self, i: usize) -> f64 { self.dihedral[i] }

    /// Sets the distance used to define `i`-th atom
    pub fn set_r(&mut self, i: usize, new_r: f64)  { self.dihedral[i] = new_r; }

    /// Sets the planar angle used to define `i`-th atom
    pub fn set_planar(&mut self, i: usize, new_planar: f64) { self.planar[i] = new_planar; }

    /// Sets the dihedral angle used to define `i`-th atom
    pub fn set_dihedral(&mut self, i: usize, new_dihedral: f64) { self.dihedral[i] = new_dihedral; }

    /// Adds an atom to this [`KinematicAtomTree`](KinematicAtomTree)
    pub fn add_atom(&mut self, atom: &InternalAtomDefinition, index: usize) {
        match self.reference_atoms.len() {
            0 => { self.reference_atoms.push([0, 0, 0, index]) }
            1 => { self.reference_atoms.push([0, 1, 0, index]) }
            2 => { self.reference_atoms.push([0, 1, 2, index]) }
            _ => {
                let a_idx = self.last_atom_by_name.get(&atom.a_name).cloned().unwrap_or(0usize);
                let b_idx = self.last_atom_by_name.get(&atom.b_name).cloned().unwrap_or(0usize);
                let c_idx = self.last_atom_by_name.get(&atom.c_name).cloned().unwrap_or(0usize);
                self.reference_atoms.push([a_idx, b_idx, c_idx, index]);
            }
        }
        self.last_atom_by_name.insert(atom.name.clone(), index);
        self.r.push(atom.r);
        self.planar.push(atom.planar);
        self.dihedral.push(atom.dihedral);
        self.names.insert(index, atom.name.clone());
    }

    /// Restores Cartesian coordinates of all the atoms of this tree
    pub fn restore_atoms(&self) -> Vec<Vec3> {
        let mut result = vec![Vec3::from_float(0.0); self.len()];
        restore_branched_chain(&self.r, &self.planar, &self.dihedral, &self.reference_atoms, &mut result);

        return result;
    }
}