use crate::coordinates::{Coordinates};
use crate::coordinates_view::{CoordinatesView};
use bioshell_numerical::Vec3;
use bioshell_sim::System;
use crate::trait_nb_list_rules::NbListRules;


impl Clone for Box<dyn NbListRules> {
    fn clone(&self) -> Box<dyn NbListRules> {
        self.box_clone()
    }
}

#[derive(Clone)]
/// Implements interaction rules for a monoatomic fluid simulation, such as argon gas.
///
/// According to these simple rules all atoms always interact which each other, so all pairs
/// are included in NBL hashing.
///
pub struct ArgonRules;

impl NbListRules for ArgonRules {
    /// Always returns ``false``, because all atoms interact.
    /// All arguments of this method are therefore ignored
    fn if_atom_excluded(&self, _coordinates: &Coordinates, _i_atom: usize) -> bool {
        false
    }

    /// Returns ``false``, because all atom pairs are included in pairwise interactions.
    ///
    /// The only exception is when ``i_atom == j_atom``: this will return false because an atom
    /// can't interact with itself
    fn if_pair_excluded(&self, _coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool {
        i_atom == j_atom
    }

    fn box_clone(&self) -> Box<dyn NbListRules> {
        Box::new((*self).clone())
    }
}

/// Implements simple polymer style interaction rules.
///
/// According to these rules:
///    - every atom is included
///    - direct neighbors do not interact as it is assumed they are connected with a pseudo-bond,
///      e.g. a harmonic spring. More precisely, ``if_pair_excluded(i, j)`` returns ``true`` if
///      ``|i-j| < 2`` and the two atoms belong to the same chain
///
#[derive(Clone)]
pub struct PolymerRules;

impl NbListRules for PolymerRules {
    /// Always returns ``false``, because all atoms of a polymer interact
    /// All arguments of this method are therefore ignored
    fn if_atom_excluded(&self, _coordinates: &Coordinates, _i_atom: usize) -> bool {
        false
    }

    /// Excludes direct neighbors in a chain.
    ///
    /// I.e. returns ``true`` if ``|i-j| < 2`` and the two atoms belong to the same chain;
    /// this prevents contact (non-bonded) interaction between atoms (beads) that are connected
    /// with bonds (e.g. harmonic springs). For any other pair of atoms this method returns `true`
    fn if_pair_excluded(&self, coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool {
        // ---------- atoms of different chains always interact
        if coordinates[i_atom].chain_id != coordinates[j_atom].chain_id {
            return false;
        }

        // ---------- exclude direct neighbors in a chain
        if i_atom > j_atom {
            return i_atom - j_atom < 2;
        } else {
            return j_atom - i_atom < 2;
        }
    }

    fn box_clone(&self) -> Box<dyn NbListRules> {
        Box::new((*self).clone())
    }
}

macro_rules! insert_nb_pair {
    ($i:expr, $j:expr, $system:expr, $self:expr) => {
        if !$self.nb_rules.if_pair_excluded(&$system, $i, $j) {
            if $system.get_closest_distance_square($j, $i) < $self.total_cutoff_sq {
                $self.nb_lists[$j].push($i);
                $self.nb_lists[$i].push($j);
            }
        }
    };
}

pub struct NeighborList {
    cutoff: f64,                    // distance cutoff for this list
    buffer_width: f64,              // safety buffer thickness
    total_cutoff_sq: f64,           // total cutoff = cutoff + safety buffer
    max_moved_sq: f64,              // how far atom can travel to trigger update
    nb_rules: Box<dyn NbListRules>, // encodes rules for this NBL
    recent_pos: Vec<Vec3>, // positions from the previous update() method call to compute displacement vector
    travelled: Vec<Vec3>, // vector of accumulative displacement of i-th atom since the last NBL update
    nb_lists: Vec<Vec<usize>>, // the NBL itself
}



impl NeighborList {
    pub fn new(cutoff: f64, buffer_thickness: f64, nb_rules: Box<dyn NbListRules>) -> NeighborList {
        let neighbors: Vec<Vec<usize>> = Vec::new();
        let total = buffer_thickness + cutoff;
        let max_moved = buffer_thickness / 2.0;
        NeighborList {
            cutoff,
            buffer_width: buffer_thickness,
            total_cutoff_sq: total * total,
            max_moved_sq: max_moved * max_moved,
            nb_rules,
            recent_pos: Vec::new(),
            travelled: Vec::new(),
            nb_lists: neighbors,
        }
    }

    pub fn get_recent_pos(&self)->&Vec<Vec3>{&self.recent_pos}
    pub fn get_travelled(&self)->&Vec<Vec3>{&self.travelled}
    pub fn get_neighbors(&self) -> &Vec<Vec<usize>> {
        &self.nb_lists
    }
    /// Provide a read-only access to neighbors of a given atom
    pub fn get_neighbor(&self, pos: usize) -> &Vec<usize> {
        &self.nb_lists[pos]
    }
    /// Provides the interaction cutoff radius
    pub fn get_cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Modifies the interaction cutoff radius
    pub fn set_cutoff(&mut self, d0: f64) {
        self.cutoff = d0;
        let total = self.buffer_width + self.cutoff;
        self.total_cutoff_sq = total * total;
    }

    /// Provides the width of the buffer zone
    pub fn get_buffer_thickness(&self) -> f64 {
        self.buffer_width
    }

    /// Modifies the width of the buffer zone
    pub fn set_buffer_thickness(&mut self, width: f64) {
        self.buffer_width = width;
        let total = self.buffer_width + self.cutoff;
        self.total_cutoff_sq = total * total;
        let max_moved = width / 2.0;
        self.max_moved_sq = max_moved * max_moved;
    }

    /// Sets new rules that define which atoms and atom pairs can be neighbors.
    /// Note: remember to update this neighbor list after this call!
    pub fn set_rules(&mut self, nb_rules: Box<dyn NbListRules>) {
        self.nb_rules = nb_rules;
    }

    /// Updates the list of neighbors after a given atom was moved.
    pub fn update(&mut self, system: &Coordinates, atom_index: usize) {
        // ---------- Is the atom relevant for this list?
        if self.nb_rules.if_atom_excluded(system, atom_index) {
            return;
        }

        // ---------- Extend the data structure if needed
        self.extend(&system);

        // ---------- accumulate the displacement
        let x = system.get_x(atom_index);
        let y = system.get_y(atom_index);
        let z = system.get_z(atom_index);
        self.travelled[atom_index].x += x - self.recent_pos[atom_index].x;
        self.travelled[atom_index].y += y - self.recent_pos[atom_index].y;
        self.travelled[atom_index].z += z - self.recent_pos[atom_index].z;

        // --- record position after move even if the list is not to be updated
        self.recent_pos[atom_index].x = x;
        self.recent_pos[atom_index].y = y;
        self.recent_pos[atom_index].z = z;

        // --- check if it moved far enough; if not - skip it
        if self.travelled[atom_index].length_squared() < self.max_moved_sq {
            return;
        }

        // --- reset the displacement because we are updating the NBL for that atom
        self.travelled[atom_index].x = 0.0;
        self.travelled[atom_index].y = 0.0;
        self.travelled[atom_index].z = 0.0;

        // --- First clear information about which_atom in its partners
        for i_nb in 0..self.nb_lists[atom_index].len() {
            // --- i_nb is the index of a neighbor
            let j: usize = self.nb_lists[atom_index][i_nb]; // --- for each j neighbor of i, find i on its list and remove
            if let Some(index) = self.nb_lists[j]
                .iter()
                .position(|value| *value == atom_index)
            {
                self.nb_lists[j].swap_remove(index);
            } else {
                println!("INCONSISTENT NBL!")
            }
        }
        // --- ... then clear neighbors of which_atom
        self.nb_lists[atom_index].clear();

        // --- detect neighbors
        for i in 0..atom_index {
            // --- check distance, update list of neighbors
            insert_nb_pair!(atom_index, i, system, self);
        }
        for i in atom_index + 1..system.get_size() {
            insert_nb_pair!(atom_index, i, system, self);
        }
    }

    /// Creates list of neighbors for each atom; the previous contents is wiped out
    pub fn update_all(&mut self, system: &Coordinates) {
        // --- extend the list if needed
        self.extend(&system);
        // --- clear all neighbors
        for vi in 0..self.nb_lists.len() {
            self.nb_lists[vi].clear()
        }
        // --- copy coordinates
        for i in 0..system.get_size() {
            if self.nb_rules.if_atom_excluded(system, i) {
                continue;
            }
            self.recent_pos[i].x = system.get_x(i);
            self.recent_pos[i].y = system.get_y(i);
            self.recent_pos[i].z = system.get_z(i);
            self.travelled[i].x = 0.0;
            self.travelled[i].y = 0.0;
            self.travelled[i].z = 0.0;
            for j in 0..i {
                // --- exclude excluded pairs, insert the relevant ones - now all in a single macro
                insert_nb_pair!(j, i, system, self);
            }
        }
    }

    /// Updates the list of neighbors after a given atom was moved.
    pub fn update_for_view(&mut self, system_view: CoordinatesView<'_>, pos: usize) {
        self.update(system_view.points, pos);
    }

    /// Creates list of neighbors for each atom; the previous contents is wiped out
    pub fn update_all_for_view(&mut self, system_view: CoordinatesView<'_>) {
        self.update_all(system_view.points);
    }



    /// Add inner vectors to this non-bonded list so it has
    /// at least as many rows as the number of atoms in the
    /// given structure
    fn extend(&mut self, system: &Coordinates) {
        if system.get_size() > self.nb_lists.len() {
            for _ in self.nb_lists.len()..system.get_size() {
                self.nb_lists.push(Vec::new());
                self.recent_pos.push(Vec3::new(0.0, 0.0, 0.0));
                self.travelled.push(Vec3::new(0.0, 0.0, 0.0));
            }
        }
    }
}

impl Clone for NeighborList {
    fn clone(&self) -> NeighborList {
        // ---------- Create a deep copy of NBL rules object
        let rules_copy: Box<dyn NbListRules> = self.nb_rules.clone();
        // ---------- Create a new non-bonded list
        let mut nbl: NeighborList = NeighborList::new(self.cutoff, self.buffer_width, rules_copy);
        // ---------- Copy the content of self list to the new copy
        for i in 0..self.nb_lists.len() {
            nbl.nb_lists.push(self.nb_lists[i].clone());
        }
        // ---------- Copy the recent coordinates
        nbl.recent_pos = self.recent_pos.clone();
        // ---------- Copy the displacement vector
        nbl.travelled = self.travelled.clone();

        return nbl;
    }
}
