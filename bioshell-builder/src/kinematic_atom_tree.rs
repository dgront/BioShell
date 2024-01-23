use std::collections::HashSet;
use std::ops::Range;
use bioshell_pdb::calc::Vec3;
use bioshell_pdb::PdbAtom;
use crate::{BuilderError, InternalAtomDefinition, RelaviveResidueLocator};
use crate::nerf::{restore_branched_chain_in_order};
use crate::BuilderError::{DefinedAtomNotFound, DihedralAngleNotFound, ResidueNotDefined};

struct ResidueAtomsDefinition {
    atoms_def: Vec<InternalAtomDefinition>,
}

impl ResidueAtomsDefinition {

    /// Returns the 3-letter code name of this residue.
    pub fn res_name(&self) -> &String {
        &self.atoms_def[0].res_name
    }

    /// Provides the definition of an atom given by its index
    pub fn atom_definition(&self, atom_index: usize) -> &InternalAtomDefinition { &self.atoms_def[atom_index] }

    /// Index of an atom given its name, starts from 0.
    pub fn atom_index(&self, atom_name: &str) -> Option<usize> {
        let mut i= 0;
        for a in &self.atoms_def {
            if a.name == atom_name { return Some(i); }
            i += 1;
        }
        return None;
    }

    /// Appends a new atom definition to this residue
    pub fn add_atom(&mut self, atom_def: &InternalAtomDefinition) { self.atoms_def.push(atom_def.clone())}

    /// Creates a new, empty residue of a given name.
    pub fn new(name: &str) -> ResidueAtomsDefinition {
        ResidueAtomsDefinition{ atoms_def: vec![] }
    }

    /// Number of atoms defined for this residue
    pub fn len(&self) -> usize { self.atoms_def.len() }
}

/// Computes Cartesian coordinates of atoms from theirs internal definition
pub struct KinematicAtomTree {
    defined_residues: Vec<ResidueAtomsDefinition>,
    r: Vec<f64>,
    planar: Vec<f64>,
    dihedral: Vec<f64>,
    names: Vec<String>,
    elements: Vec<String>,
    dihedral_names: Vec<String>,
    reference_atoms: Vec<[usize;4]>,
    building_order: Vec<usize>,
    residue_atoms: Vec<Range<usize>>,
    is_compiled: bool
}

impl KinematicAtomTree {

    /// Creates an empty [`KinematicAtomTree`](KinematicAtomTree).
    pub fn new() -> KinematicAtomTree {
        KinematicAtomTree {
            defined_residues: vec![],
            r: vec![],
            planar: vec![],
            dihedral: vec![],
            names: Default::default(),
            elements: vec![],
            dihedral_names: vec![],
            reference_atoms: vec![],
            building_order: vec![],
            residue_atoms: vec![],
            is_compiled: false
        }
    }

    /// Counts atoms of this polymer
    pub fn count_atoms(&mut self) -> usize {
        if !self.is_compiled { self.build_internal_data(); }
        self.names.len()
    }

    /// Counts residues of this polymer
    pub fn count_residues(&mut self) -> usize { self.defined_residues.len() }

    pub fn atoms_for_residue(&mut self, residue_index: usize) -> &Range<usize> {
        if !self.is_compiled { self.build_internal_data(); }
        &self.residue_atoms[residue_index]
    }

    pub fn residue_for_atom(&self, atom_index: usize) -> usize {
        for i_res in 0..self.residue_atoms.len() {
            if atom_index >= self.residue_atoms[i_res].start && atom_index < self.residue_atoms[i_res].end {
                return i_res;
            }
        }
        panic!("Atom index out of bounds");
    }

    /// Provides a name for an atom given its index
    pub fn atom_name(&mut self, atom_index: usize) -> &String {
        if !self.is_compiled { self.build_internal_data(); }
        &self.names[atom_index]
    }

    /// Provides a chemical element for an atom given its index
    pub fn atom_element(&mut self, atom_index: usize) -> &String {
        if !self.is_compiled { self.build_internal_data(); }
        &self.elements[atom_index]
    }

    ///  Provides a name for a residue given its index.
    ///
    /// Note that the index of a residue for a given atom may be found by calling
    /// [`residue_for_atom()`](KinematicAtomTree::residue_for_atom) method.
    pub fn residue_name(&mut self, residue_index: usize)  -> &String {
        return self.defined_residues[residue_index].res_name();
    }

    /// Provides a dihedral angle value given the name of that angle.
    ///
    /// The dihedral angle name must match the name used by [`InternalAtomDefinition`](InternalAtomDefinition) of the respective atom.
    /// Returns None when the requested dihedral angle hasn't been defined for  given re
    pub fn named_dihedral(&mut self, residue_index: usize, dihedral_name: &str) -> Option<f64> {
        if !self.is_compiled { self.build_internal_data(); }
        for i_atom in self.residue_atoms[residue_index].start..self.residue_atoms[residue_index].end {
            if self.dihedral_names[i_atom] == dihedral_name {
                return Some(self.dihedral[i_atom]);
            }
        }
        return None;
    }

    /// Sets a new value for a dihedral angle given the name of that angle.
    ///
    /// The dihedral angle name must match the name used by [`InternalAtomDefinition`](InternalAtomDefinition) of the respective atom.
    pub fn set_named_dihedral(&mut self, residue_index: usize, dihedral_name: &str, value: f64) -> Result<(), BuilderError> {
        if !self.is_compiled { self.build_internal_data()?; }
        for i_atom in self.residue_atoms[residue_index].start..self.residue_atoms[residue_index].end {
            if self.dihedral_names[i_atom] == dihedral_name {
                self.dihedral[i_atom] = value;
                return Ok(());
            }
        }
        return Err(DihedralAngleNotFound { residue_index, dihedral_name: dihedral_name.to_string() })
    }

    pub fn add_atom(&mut self, atom: &InternalAtomDefinition, residue_index: usize) {
        while self.defined_residues.len() <= residue_index {
            self.defined_residues.push(ResidueAtomsDefinition::new("") );
        }
        self.defined_residues[residue_index].add_atom(atom);
        self.is_compiled = false;
    }

    pub fn add_residue(&mut self, residue_definition: &Vec<InternalAtomDefinition>) {
        let residue_index = self.defined_residues.len();
        for atom in residue_definition {
            self.add_atom(atom, residue_index);
        }
        self.is_compiled = false;
    }

    pub fn patch_residue(&mut self, residue_index: usize, residue_definition: &Vec<InternalAtomDefinition>) -> Result<(), BuilderError> {
        if residue_index >= self.defined_residues.len() {
            return Err(ResidueNotDefined { residue_index: residue_index })
        }

        for atom_def in residue_definition {
            if let Some(atom_index) = self.defined_residues[residue_index].atom_index(&atom_def.name) {
                self.defined_residues[residue_index].atoms_def[atom_index] = atom_def.clone();
            } else {
                self.defined_residues[residue_index].add_atom(atom_def);
            }
        }
        self.is_compiled = false;
        return Ok(());
    }

    fn atom_index(&self, current_residue_index: usize, relative_residue_location: &RelaviveResidueLocator,
                  atom_name: &str) -> Result<usize, BuilderError> {
        let the_residue = match relative_residue_location {
            RelaviveResidueLocator::Previous => { current_residue_index - 1 }
            RelaviveResidueLocator::This => { current_residue_index }
            RelaviveResidueLocator::Next => { current_residue_index + 1 }
        };
        if the_residue >= self.defined_residues.len() {
            return Err(ResidueNotDefined { residue_index: the_residue })
        }
        let atom_local_index = self.defined_residues[the_residue].atom_index(atom_name);
        if atom_local_index.is_none() {
            return Err(DefinedAtomNotFound { residue_index: the_residue, atom_name: atom_name.to_string() })
        }
        let mut atoms_total = 0;
        if the_residue > 0 {
            for ires in 0..the_residue {
                atoms_total += self.defined_residues[ires].len();
            }
        }

        return Ok(atom_local_index.unwrap() + atoms_total);
    }

    fn setup_residue_ranges(&mut self) {
        self.residue_atoms = vec![0..0; self.count_residues()];
        self.residue_atoms[0] = 0..self.defined_residues[0].len();
        for i_res in 1..self.defined_residues.len() {
            let n = self.defined_residues[i_res].len();
            self.residue_atoms[i_res] = self.residue_atoms[i_res-1].end..self.residue_atoms[i_res-1].end+n;
        }
    }

    /// Run setup_residue_ranges() method before calling this one!
    fn build_internal_data(&mut self) -> Result<(), BuilderError> {
        let n_atoms: usize = self.defined_residues.iter().map(|v| v.len()).sum();
        self.r.resize(n_atoms, 0.0);
        self.planar.resize(n_atoms, 0.0);
        self.names.resize(n_atoms, Default::default());
        self.elements.resize(n_atoms, Default::default());
        self.dihedral_names.resize(n_atoms, Default::default());
        self.dihedral.resize(n_atoms, 0.0);
        self.building_order.clear();
        self.reference_atoms.resize(n_atoms, [0, 0, 0, 0]);
        let mut i_atom: usize = 0;
        let mut i_residue: usize = 0;

        for residue in &self.defined_residues {

            for i_atom_def in 0..residue.len() {
                let atom = residue.atom_definition(i_atom_def);
                self.r[i_atom] = atom.r;
                self.planar[i_atom] = atom.planar;
                self.dihedral[i_atom] = atom.dihedral;
                self.names[i_atom] = atom.name.clone();
                self.dihedral_names[i_atom] = atom.dihedral_name.clone();
                self.elements[i_atom] = atom.element.clone();
                match i_atom {
                    0 => self.reference_atoms[i_atom] = [0, 0, 0, 0],
                    1 => self.reference_atoms[i_atom] = [0, 0, 0, 1],
                    2 => self.reference_atoms[i_atom] = [0, 0, 1, 2],
                    _ => {
                        let a = self.atom_index(i_residue, &atom.a_residue, &atom.a_name);
                        let b = self.atom_index(i_residue, &atom.b_residue, &atom.b_name);
                        let c = self.atom_index(i_residue, &atom.c_residue, &atom.c_name);
                        let d = self.atom_index(i_residue, &atom.d_residue, &atom.name);
                        self.reference_atoms[i_atom] = [a?, b?, c?, d?]
                    }
                }
                i_atom += 1;
            }
            i_residue += 1;
        }
        // ---------- Change the order in which the atoms are reconstructed
        // ---------- Atom definition may be based on a future atom
        let mut max_atom_built = 0;
        let mut waiting_list: HashSet<(usize, usize)> = Default::default(); // --- (max_atom_index_required, atom_index_built)
        let mut waiting_inserted: HashSet<(usize, usize)> = Default::default();
        for i in 0..self.reference_atoms.len() {
            // ---------- First check if we can insert any atom from the waiting list
            for el in &waiting_list {
                if el.0 <= max_atom_built {
                    self.building_order.push(el.1);
                    max_atom_built = el.1;
                    waiting_inserted.insert((el.0, el.1));
                }
            }
            for el in &waiting_inserted { waiting_list.remove(&(el.0, el.1));}
            waiting_inserted.clear();
            // ---------- check maximum index of an atom used in the definition of the i-th atom
            let max_ref_atom_required = *self.reference_atoms[i].iter().max().unwrap();
            if max_ref_atom_required <= self.reference_atoms[i][3] {
                self.building_order.push(self.reference_atoms[i][3]);
                max_atom_built = max_atom_built.max(self.reference_atoms[i][3]);
            }
            else {
                waiting_list.insert((max_ref_atom_required, self.reference_atoms[i][3]));
            }
        }
        for el in &waiting_list {
            self.building_order.push(el.1);
        }

        self.setup_residue_ranges();
        self.is_compiled = true;
        return Ok(());
    }


    /// Restores Cartesian coordinates of all the atoms of this tree
    pub fn build_coordinates(&mut self) -> Result<Vec<Vec3>, BuilderError> {
        let result = self.build_internal_data();
        if result.is_err() { return Err(result.unwrap_err())}
        let mut result = vec![Vec3::from_float(0.0); self.reference_atoms.len()];
        restore_branched_chain_in_order(&self.r, &self.planar, &self.dihedral, &self.reference_atoms,
                                        &self.building_order, &mut result);

        return Ok(result);
    }

    /// Restores Cartesian coordinates of all the atoms of this tree
    pub fn build_atoms(&mut self, chain_id: &str) -> Result<Vec<PdbAtom>, BuilderError> {

        if !self.is_compiled {let _ = self.build_internal_data()?;}

        let mut vectors = vec![Vec3::from_float(0.0); self.reference_atoms.len()];
        restore_branched_chain_in_order(&self.r, &self.planar, &self.dihedral, &self.reference_atoms,
                                        &self.building_order, &mut vectors);

        let mut atoms = vec![PdbAtom::new(); self.reference_atoms.len()];
        for i in 0..atoms.len() {
            let res_seq = self.residue_for_atom(i);
            atoms[i].serial = (i + 1) as i32;
            atoms[i].res_seq = res_seq as i32 + 1;
            atoms[i].pos = vectors[i];
            atoms[i].res_name = self.residue_name(res_seq).clone();
            atoms[i].chain_id = chain_id.to_string();
            atoms[i].name = self.atom_name(i).clone();
            atoms[i].element = Some(self.elements[i].clone());
        }

        return Ok(atoms);
    }
}

