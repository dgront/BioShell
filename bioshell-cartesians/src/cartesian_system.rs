use crate::{NbList, NbListRules};
use crate::{Coordinates, CoordinatesView};
use bioshell_sim::System;

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

    pub fn set(&mut self, i:usize, x: f64, y: f64, z: f64) {
        self.coordinates.set(i, x, y, z);
        self.neighbor_list.update(&self.coordinates, i);
    }

    pub fn add(&mut self, i:usize, x: f64, y: f64, z: f64) {
        self.coordinates.add(i, x, y, z);
        self.neighbor_list.update(&self.coordinates, i);
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
    /// Note: remember to update this neighbor list after this call!
    pub fn set_rules(&mut self, nb_rules: Box<dyn NbListRules>) {
        self.neighbor_list.set_rules(nb_rules);
        let stls_v = CoordinatesView { points: &self.coordinates, };
        self.neighbor_list.update_all_for_view(stls_v);
    }
}

impl System for CartesianSystem {
    /// Returns the number of atoms of this system
    fn size(&self) -> usize { self.coordinates.size() }


    fn copy_from(&mut self, i:usize, rhs: &CartesianSystem) {
        self.coordinates.copy(i,&rhs.coordinates());
        let stls_v = CoordinatesView { points: &self.coordinates, };
        self.neighbor_list.update_for_view(stls_v, i);
    }
}