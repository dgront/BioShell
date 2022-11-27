use bioshell_numerical::{Vec3};
use bioshell_core::structure::{Coordinates, CoordinatesView};


/// Rules that define which atoms and which atom pairs will be excluded from hashing
pub trait NbListRules {
    fn if_atom_excluded(&self, coordinates: &Coordinates, i_atom:usize) -> bool;
    fn if_pair_excluded(&self, coordinates: &Coordinates, i_atom:usize, j_atom:usize) -> bool;
    /// Each NbListRules must provide a way to clone its boxed instance
    fn box_clone(&self) -> Box<dyn NbListRules>;
}

impl Clone for Box<dyn NbListRules>  {
    fn clone(&self) -> Box<dyn NbListRules> { self.box_clone() }
}

#[derive(Clone)]
pub struct ArgonRules;

impl NbListRules for ArgonRules {

    fn if_atom_excluded(&self, _coordinates: &Coordinates, _i_atom: usize) -> bool { false }

    fn if_pair_excluded(&self, _coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool {
        i_atom == j_atom
    }

    fn box_clone(&self) -> Box<dyn NbListRules> {
        Box::new((*self).clone())
    }
}

#[derive(Clone)]
pub struct PolymerRules;

impl NbListRules for PolymerRules {

    fn if_atom_excluded(&self, _coordinates: &Coordinates, _i_atom: usize) -> bool { false }

    fn if_pair_excluded(&self, _coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool {
        if i_atom >  j_atom {
            return i_atom - j_atom < 2;
        } else {
            return j_atom - i_atom < 2;
        }
    }

    fn box_clone(&self) -> Box<dyn NbListRules> {
        Box::new((*self).clone())
    }
}

pub struct NbList {
    cutoff: f64,                        // distance cutoff for this list
    buffer_width: f64,                  // safety buffer thickness
    total_cutoff_sq: f64,               // total cutoff = cutoff + safety buffer
    max_moved_sq: f64,                  // how far atom can travel to trigger update
    nb_rules: Box<dyn NbListRules>,     // encodes rules for this NBL
    recent_pos: Vec<Vec3>,              // positions from last NBL update
    nb_lists: Vec<Vec<usize>>,          // the NBL itself
}

macro_rules! insert_nb_pair {
    ($i:expr, $j:expr, $system:expr, $self:expr) => {
        if ! $self.nb_rules.if_pair_excluded(&$system, $i, $j) {
            if $system.closest_distance_square($j,$i) < $self.total_cutoff_sq {
                $self.nb_lists[$j].push($i);
                $self.nb_lists[$i].push($j);
            }
        }
    }
}

impl NbList {

    pub fn new(cutoff:f64, buffer_thickness: f64, nb_rules: Box<dyn NbListRules>) -> NbList {
        let neighbors:Vec<Vec<usize>> = Vec::new();
        let total = buffer_thickness + cutoff;
        let max_moved = buffer_thickness / 2.0;
        NbList { cutoff, buffer_width: buffer_thickness, total_cutoff_sq: total*total,
            max_moved_sq:max_moved*max_moved, nb_rules, recent_pos: Vec::new(), nb_lists: neighbors }
    }

    /// Provides the interaction cutoff radius
    pub fn cutoff(&self) -> f64 { self.cutoff }

    /// Modifies the interaction cutoff radius
    pub fn set_cutoff(&mut self, d0: f64) {
        self.cutoff = d0;
        let total = self.buffer_width + self.cutoff;
        self.total_cutoff_sq = total*total;
    }

    /// Provides the width of the buffer zone
    pub fn buffer_width(&self) -> f64 { self.buffer_width }

    /// Modifies the width of the buffer zone
    pub fn set_buffer_width(&mut self, width: f64) {
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
    pub fn update(&mut self, system: &Coordinates, which_atom: usize) {

        // --- Is the atom relevant for this list?
        if self.nb_rules.if_atom_excluded(system, which_atom) { return;}
        // --- extend the list if needed
        // self.extend(&system);
        // --- check if it moved far enough; if not - skip it
        if system.closest_distance_square_to_vec(which_atom,&self.recent_pos[which_atom]) < self.max_moved_sq {return;}
        // --- copy coordinates of the atom we are updating
        self.recent_pos[which_atom].x = system.x(which_atom);
        self.recent_pos[which_atom].y = system.y(which_atom);
        self.recent_pos[which_atom].z = system.z(which_atom);

        // --- First clear information about which_atom in its partners
        for i_nb in  0..self.nb_lists[which_atom].len() {       // --- i_nb is the index of a neighbor
            let j:usize = self.nb_lists[which_atom][i_nb];             // --- for each j neighbor of i, find i on its list and remove
            if let Some(index) = self.nb_lists[j].iter().position(|value| *value == which_atom) {
                self.nb_lists[j].swap_remove(index);
            } else { println!("INCONSISTENT NBL!") }
        }
        // --- ... then clear neighbors of which_atom
        self.nb_lists[which_atom].clear();

        // --- detect neighbors
        for i in 0..which_atom { // --- check distance, update list of neighbors
            insert_nb_pair!(which_atom, i, system, self);
        }
        for i in which_atom+1..system.size() {
            insert_nb_pair!(which_atom, i, system, self);
        }
    }

    /// Creates list of neighbors for each atom; the previous contents is wiped out
    pub fn update_all(&mut self, system: &Coordinates) {
        // --- extend the list if needed
        self.extend(&system);
        // --- clear all neighbors
        for vi in 0..self.nb_lists.len() { self.nb_lists[vi].clear() }
        // --- copy coordinates
        for  i in 0..system.size() {
            if self.nb_rules.if_atom_excluded(system, i) { continue }
            self.recent_pos[i].x = system.x(i);
            self.recent_pos[i].y = system.y(i);
            self.recent_pos[i].z = system.z(i);
            for j in 0..i {
                // --- exclude excluded pairs
                // if self.nb_rules.if_pair_excluded(i, j) { continue; }
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

    /// Provide a read-only access to neighbors of a given atom
    pub fn neighbors(&self, pos:usize) -> &Vec<usize> { &self.nb_lists[pos] }

    /// Add inner vectors to this non-bonded list so it has at least as many rows as the number of atoms in the given structure
    fn extend(&mut self, system: &Coordinates) {
        if system.size() > self.nb_lists.len() {
            for _ in self.nb_lists.len()..system.size() {
                self.nb_lists.push(Vec::new());
                self.recent_pos.push(Vec3::new(0.0,0.0,0.0));
            }
        }
    }
}


impl Clone for NbList {
    fn clone(&self) -> NbList {
        // ---------- Create a deep copy of NBL rules object
        let rules_copy: Box<dyn NbListRules> = self.nb_rules.clone();
        // ---------- Create a new non-bonded list
        let mut nbl: NbList = NbList::new(self.cutoff, self.buffer_width, rules_copy);
        // ---------- Copy the content of self list to the new copy
        for i in 0..self.nb_lists.len() {
            nbl.nb_lists.push(self.nb_lists[i].clone());
        }
        // ---------- Copy the recent coordinates
        nbl.recent_pos = self.recent_pos.clone();

        return nbl;
    }
}