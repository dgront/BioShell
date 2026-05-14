use std::collections::HashMap;
use petgraph::graph::{Graph, NodeIndex};
use petgraph::prelude::Undirected;
use petgraph::visit::EdgeRef;

use crate::{Atom, BondType, ChemErrors};


/// A molecule is a graph of atoms connected with bonds.
///
/// Internally, atoms are stored as graph vertices and bonds are stored as
/// graph edges with `BondType` as edge data.
///
/// # Examples
/// ```
/// # use bioshell_chem::{Atom, BondType, ChemErrors, Molecule};
///
/// // Creates a benzene molecule with 6 carbon atoms and aromatic bonds.
/// fn benzene() -> Result<Molecule, ChemErrors> {
///     let mut mol = Molecule::new("benzene"); // Empty molecule
///     for i in 0..6 {                         // Add 6 carbon atoms with indices from 0 to 5
///         mol.add_atom(Atom::new(i, 6, 0))?;
///     }
///
///     for i in 0..6 {
///         mol.bind_atoms(i, (i + 1) % 6, BondType::Aromatic)?;
///     }
///     Ok(mol)
/// }
///
/// # fn main() -> Result<(), ChemErrors> {
/// let mol = benzene()?;
/// assert_eq!(mol.molecule_name, "benzene");
/// assert_eq!(mol.count_atoms(), 6);
/// assert_eq!(mol.count_bonds(), 6);
/// # Ok(())
/// # }
/// ```
pub struct Molecule {
    /// The name of this molecule.
    ///
    /// Molecule name is required by MOL2 file format from TRIPOS.
    pub molecule_name: String,

    graph: Graph<Atom, BondType, Undirected>,
    vertex_flags: Vec<u8>,
    atom_id_to_node: HashMap<usize, NodeIndex>,
}


impl Molecule {

    /// Creates a new molecule with the given name.
    ///
    /// # Examples
    /// ```
    /// use bioshell_chem::Molecule;
    /// let mol = Molecule::new("benzene");
    /// assert_eq!(mol.molecule_name, "benzene");
    /// assert_eq!(mol.count_atoms(), 0);
    /// assert_eq!(mol.count_bonds(), 0);
    /// ```
    pub fn new(name: &str) -> Self {
        Self {
            molecule_name: name.to_string(),
            graph: Graph::new_undirected(),
            vertex_flags: Vec::new(),
            atom_id_to_node: HashMap::default(),
        }
    }

    /// Returns the number of atoms in this molecule.
    pub fn count_atoms(&self) -> usize {
        self.graph.node_count()
    }

    /// Adds a new atom to this molecule.
    ///
    /// The atom is identified by its internal index, which can be used to preserve the original atom numbering
    pub fn add_atom(&mut self, atom: Atom) -> Result<(), ChemErrors>  {
        let id = atom.index();

        if self.atom_id_to_node.contains_key(&id) {
            return Err(ChemErrors::DuplicateAtomId(id));
        }

        let node_index = self.graph.add_node(atom);
        self.atom_id_to_node.insert(id, node_index);
        self.vertex_flags.push(0);

        Ok(())
    }

    /// Returns true if the given atom index belongs to this molecule.
    pub fn contains_atom(&self, atom_idx: usize) -> bool {
        self.atom_id_to_node.contains_key(&atom_idx)
    }

    /// Returns an atom referred by an index.
    pub fn get_atom(&self, atom_idx: usize) -> Option<&Atom> {
        if let idx = self.atom_id_to_node[&atom_idx] {
            self.graph.node_weight(idx)
        } else {
            None
        }
    }

    /// Returns a mutable atom referred by an index.
    // pub fn get_atom_mut(&mut self, atom_idx: usize) -> Option<&mut Atom> {
    //     if let idx = self.atom_id_to_node[&atom_idx] {
    //         self.graph.node_weights_mut(idx)
    //     } else {
    //         None
    //     }
    // }

    /// Returns an iterator over atoms in this molecule.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.graph.node_weights()
    }

    /// Returns a mutable iterator over atoms in this molecule.
    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut Atom> {
        self.graph.node_weights_mut()
    }

    /// Binds two atoms with a bond.
    ///
    /// If the bond already exists, its `BondType` is replaced.
    pub fn bind_atoms(&mut self, first_atom_idx: usize, second_atom_idx: usize, bond_type: BondType) -> Result<(), ChemErrors> {
        let first = self.atom_id_to_node[&first_atom_idx];
        let second = self.atom_id_to_node[&second_atom_idx];

        self.graph.update_edge(first, second, bond_type);
        Ok(())
    }

    /// Removes a bond between two given atoms.
    ///
    /// Returns true if the bond was actually removed.
    pub fn remove_bond(&mut self, first_atom_idx: usize, second_atom_idx: usize) -> Result<bool, ChemErrors> {
        let first_atom = self.atom_id_to_node[&first_atom_idx];
        let second_atom = self.atom_id_to_node[&second_atom_idx];

        if let Some(edge_index) = self.graph.find_edge(first_atom, second_atom) {
            self.graph.remove_edge(edge_index);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Returns true if the two atoms are connected with a bond.
    pub fn are_bonded(&self, first_atom_idx: usize, second_atom_idx: usize) -> Result<bool, ChemErrors> {
        let first_atom = self.atom_id_to_node[&first_atom_idx];
        let second_atom = self.atom_id_to_node[&second_atom_idx];

        Ok(self.graph.find_edge(first_atom, second_atom).is_some())
    }

    /// Returns the number of bonds in this molecule.
    pub fn count_bonds(&self) -> usize {
        self.graph.edge_count()
    }

    /// Returns the number of bonds connected to the given atom in this molecule.
    pub fn count_bonds_for_atom(&self, atom_idx: usize) -> usize {
        let atom = self.atom_id_to_node[&atom_idx];
        self.graph.neighbors(atom).count()
    }

    /// Returns a bond between two atoms.
    pub fn get_bond(&self, first_atom_idx: usize, second_atom_idx: usize) -> Result<&BondType, ChemErrors> {
        let first_atom = self.atom_id_to_node[&first_atom_idx];
        let second_atom = self.atom_id_to_node[&second_atom_idx];

        let edge_index = self
            .graph
            .find_edge(first_atom, second_atom)
            .ok_or_else(|| {
                ChemErrors::BondNotFound(first_atom.index(), second_atom.index())
            })?;

        Ok(&self.graph[edge_index])
    }

    /// Returns an iterator over bonds of this molecule.
    ///
    /// Each item is `(first_atom_id, second_atom_id, bond_type)`.
    pub fn bonds(&self) -> impl Iterator<Item = (usize, usize, &BondType)> {
        self.graph.edge_references().map(|edge| {
            let first_atom_id = self.graph[edge.source()].index();
            let second_atom_id = self.graph[edge.target()].index();

            (first_atom_id, second_atom_id, edge.weight())
        })
    }

    /// Returns atoms connected to a given atom.
    pub fn neighbors(&self, atom_idx: usize) -> Vec<&Atom> {
        let atom = self.atom_id_to_node[&atom_idx];

        self
            .graph
            .neighbors(atom)
            .map(|neighbor| &self.graph[neighbor])
            .collect()
    }

    /// Returns an iterator over atom IDs adjacent to the given atom ID.
    pub fn neighbor_indices(&self, atom_id: usize) -> impl Iterator<Item = usize> + '_ {
        let atom = self.atom_id_to_node[&atom_id];

        self.graph
            .neighbors(atom)
            .map(|neighbor| self.graph[neighbor].index())
    }

    /// Returns a flag assigned to a vertex identified by its internal index.
    ///
    /// By default, all flags are 0.
    pub fn vertex_flag(&self, atom_idx: usize) -> Result<u8, ChemErrors> {
        let atom = self.atom_id_to_node[&atom_idx];
        Ok(self.vertex_flags[atom.index()])
    }

    /// Assigns a flag to a vertex identified by its internal index.
    pub fn set_vertex_flag(&mut self, atom_idx: usize, flag: u8) -> Result<(), ChemErrors> {
        let atom = self.atom_id_to_node[&atom_idx];
        self.vertex_flags[atom.index()] = flag;
        Ok(())
    }

    // fn check_atom_index(&self, atom: AtomIndex) -> Result<(), ChemErrors> {
    //     if self.graph.node_weight(atom).is_some() {
    //         Ok(())
    //     } else {
    //         Err(ChemErrors::AtomIndexOutOfBounds(atom.index()))
    //     }
    // }
}