use bioshell_numerical::Vec3;
use crate::{NbList, NbListRules};
use crate::{Coordinates, CoordinatesView};
use bioshell_sim::{ResizableSystem, System};

#[derive(Clone)]
pub struct CartesianSystem {
    coordinates: Coordinates,
    neighbor_list: NbList,
}

impl CartesianSystem {
    pub fn new(coordinates:Coordinates, neighbor_list:NbList) -> CartesianSystem {

        let mut out = CartesianSystem{coordinates, neighbor_list};
        out.neighbor_list.update_all(&out.coordinates);

        return out;
    }

    /// Provide immutable access to coordinates of this system
    pub fn coordinates(&self) -> &Coordinates { &self.coordinates }

    /// Current volume of the simulation box
    pub fn volume(&self) -> f64 { self.coordinates.box_len().powf(3.0) }

    /// Returns the simulation box length
    #[inline(always)]
    pub fn box_len(&self) -> f64 { self.coordinates.box_len() }

    /// Changes the simulation box length which results in the change of all positions and the volume
    pub fn set_box_len(&mut self, new_box_len: f64) {

        // ---------- expansion / contraction factor
        let f = new_box_len/self.coordinates().box_len();
        // ---------- set the new box length
        self.coordinates.set_box_len(new_box_len);
        // ---------- alter atomic positions
        for i in 0..self.size() {
            let x = self.coordinates.x(i) * f;
            let y = self.coordinates.y(i) * f;
            let z = self.coordinates.z(i) * f;
            self.coordinates.set(i, x, y, z);
        }
    }

    /// Provide immutable access to the list of neighbors
    pub fn neighbor_list(&self) -> & NbList { & self.neighbor_list }

    /// Recalculate neighbors for the i-th atom of this system.
    ///
    /// Call this method after the i-th atom has been moved.
    ///
    /// # Arguments
    /// * `i` - index of an atom to be recalculated
    pub fn update_nbl(&mut self, pos: usize) {
        self.neighbor_list.update(&self.coordinates, pos);
    }

    /// Changes this system by setting new coordinates for one of its atoms.
    ///
    /// This method **does not** trigger a non-bonded list update. To recalculate neighbors for
    /// the moved atom, call `update_nbl()` method
    ///
    /// # Arguments
    /// * `i` - index of an atom to be modified
    /// * `x` - the new X coordinate value
    /// * `y` - the new X coordinate value
    /// * `z` - the new X coordinate value
    pub fn set(&mut self, i:usize, x: f64, y: f64, z: f64) {
        self.coordinates.set(i, x, y, z);
    }

    /// Assigns the residue type for the atom `i` to `t`
    pub fn set_res_type(&mut self, i:usize, t: u8) { self.coordinates[i].res_type = t; }

    /// Assigns the type of the atom `i` to `t`
    pub fn set_atom_type(&mut self, i:usize, t: u8) { self.coordinates[i].atom_type = t; }


    /// Assign each atom of this system to a chain.
    ///
    /// This method calls [`Coordinates::set_chains()`] to assign chains of this system
    pub fn set_chains(&mut self, chain_ranges: &Vec<(usize, usize)>) {
        self.coordinates.set_chains(chain_ranges);
    }

    /// Adds x, y, z to a given atom of this system
    ///
    /// This method **does not** trigger a non-bonded list update. To recalculate neighbors for
    /// the moved atom, call `update_nbl()` method
    ///
    /// # Arguments
    /// * `i` - index of an atom to be modified
    /// * `x` - the new X coordinate value
    /// * `y` - the new X coordinate value
    /// * `z` - the new X coordinate value
    pub fn add(&mut self, i:usize, x: f64, y: f64, z: f64) {
        self.coordinates.add(i, x, y, z);
    }

    /// Copies coordinates of a given atom from a given vector.
    ///
    /// This method also triggers neighbor list update at position `i`.
    /// # Arguments
    /// * `i` - index of an atom to be copied; it's the same index in both source and destination coordinates
    /// * `rhs` - the source vector to copy x, y and z from
    pub fn copy_from_vec(&mut self, i:usize, rhs: &Vec3) {
        self.coordinates.copy_from_vec(i, rhs);
        let stls_v = CoordinatesView { points: &self.coordinates, };
        self.neighbor_list.update_for_view(stls_v, i);
    }

    /// Provides the interaction cutoff radius used by the neighbor list
    pub fn cutoff(&self) -> f64 { self.neighbor_list.cutoff() }

    /// Modifies the interaction cutoff radius used by the neighbor list
    pub fn set_cutoff(&mut self, d0: f64) { self.neighbor_list.set_cutoff(d0); }

    /// Provides the width of the buffer zone used by the neighbor list
    pub fn buffer_width(&self) -> f64 { self.neighbor_list.buffer_width() }

    /// Modifies the width of the buffer zone used by the neighbor list
    pub fn set_buffer_width(&mut self, width: f64) { self.neighbor_list.set_buffer_width(width); }

    /// Sets new rules that define which atoms and atom pairs can be neighbors.
    /// This method also triggers neighbor list update.
    pub fn set_rules(&mut self, nb_rules: Box<dyn NbListRules>) {
        self.neighbor_list.set_rules(nb_rules);
        let stls_v = CoordinatesView { points: &self.coordinates, };
        self.neighbor_list.update_all_for_view(stls_v);
    }
}

impl System for CartesianSystem {
    /// Returns the number of atoms of this system
    fn size(&self) -> usize { self.coordinates.size() }

    /// Copies coordinates of a given atom from another system.
    ///
    /// This method also triggers neighbor list update at position `i`.
    /// # Arguments
    /// * `i` - index of an atom to be copied; it's the same index in both source and destination coordinates
    /// * `rhs` - the source to copy from
    fn copy_from(&mut self, i:usize, rhs: &Self) {
        self.coordinates.copy_from(i, rhs.coordinates());
        let stls_v = CoordinatesView { points: &self.coordinates, };
        self.neighbor_list.update_for_view(stls_v, i);
    }
}

impl ResizableSystem for CartesianSystem {
    /// Change the number of atoms of this system
    fn set_size(&mut self, new_size: usize) { self.coordinates.set_size(new_size); }

    /// Returns the maximum number of atoms system may have
    fn capacity(&self) -> usize { self.coordinates.capacity() }
}