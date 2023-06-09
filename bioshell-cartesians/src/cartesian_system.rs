use crate::coordinates::{Coordinates};
use crate::coordinates_view::{CoordinatesView};
use crate::nblist::{NbList};
use crate::trait_nb_list_rules::{NbListRules};
use bioshell_numerical::Vec3;
use bioshell_sim::{ResizableSystem, System};

#[derive(Clone)]
pub struct CartesianSystem {
    coords_: Coordinates,
    neighbor_list_: NbList,
}

impl CartesianSystem {
    pub fn new(coords: Coordinates, neighbor_list: NbList) -> CartesianSystem {
        let mut system = CartesianSystem {
            coords_: coords,
            neighbor_list_: neighbor_list,
        };
        system.neighbor_list_.update_all(&system.coords_);
        return system;
    }

    /// Provide immutable access to coordinates of this system
    pub fn get_coordinates(&self) -> &Coordinates {
        &self.coords_
    }

    /// Current volume of the simulation box
    pub fn get_box_volume(&self) -> f64 {
        self.coords_.get_box_len().powf(3.0)
    }

    /// Returns the simulation box length
    #[inline(always)]
    pub fn get_box_len(&self) -> f64 {
        self.coords_.get_box_len()
    }

    /// Changes the simulation box length which results in the change of all positions and the volume
    pub fn set_box_len(&mut self, new_box_len: f64) {
        // ---------- expansion / contraction factor
        let f = new_box_len / self.get_coordinates().get_box_len();
        // ---------- set the new box length
        self.coords_.set_box_len(new_box_len);
        // ---------- alter atomic positions
        for i in 0..self.get_size() {
            let x = self.coords_.get_x(i) * f;
            let y = self.coords_.get_y(i) * f;
            let z = self.coords_.get_z(i) * f;
            self.coords_.set_xyz(i, x, y, z);
        }
    }

    /// Provide immutable access to the list of neighbors
    pub fn get_neighbor_list(&self) -> &NbList {
        &self.neighbor_list_
    }

    /// Recalculate neighbors for the i-th atom of this system.
    ///
    /// Call this method after the i-th atom has been moved.
    ///
    /// # Arguments
    /// * `i` - index of an atom to be recalculated
    pub fn update_nbl(&mut self, pos: usize) {

        self.neighbor_list_.update(&self.coords_, pos);
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
    pub fn set_xyz(&mut self, i: usize, x: f64, y: f64, z: f64) {
        self.coords_.set_xyz(i, x, y, z);
    }
    pub fn set_vec3(&mut self, i: usize, vec: Vec3) {
        self.coords_.set_vec3(i, vec);
    }

    /// Assigns the residue type for the atom `i` to `t`
    pub fn set_res_type(&mut self, i: usize, t: u8) {
        self.coords_[i].res_type = t;
    }

    /// Assigns the type of the atom `i` to `t`
    pub fn set_atom_type(&mut self, i: usize, t: u8) {
        self.coords_[i].atom_type = t;
    }

    /// Assign each atom of this system to a chain.
    ///
    /// This method calls [`Coordinates::set_chains()`] to assign chains of this system
    //pub fn set_chains(&mut self, chain_collection: &Vec<Chain>) {
    //    self.chains_vec_.set_chains(chain_collection);
    //}

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
    pub fn add_xyz(&mut self, i: usize, x: f64, y: f64, z: f64) {
        self.coords_.add_xyz(i, x, y, z);
    }

    /// Copies coordinates of a given atom from a given vector.
    ///
    /// This method also triggers neighbor list update at position `i`.
    /// # Arguments
    /// * `i` - index of an atom to be copied; it's the same index in both source and destination coordinates
    /// * `rhs` - the source vector to copy x, y and z from
    pub fn copy_from_vec(&mut self, i: usize, rhs: &Vec3) {
        self.coords_.from_vec(i, rhs);
        let stls_v = CoordinatesView {
            points: &self.coords_,
        };
        self.neighbor_list_.update_for_view(stls_v, i);
    }

    /// Provides the interaction cutoff radius used by the neighbor list
    pub fn get_cutoff(&self) -> f64 {
        self.neighbor_list_.cutoff()
    }

    /// Modifies the interaction cutoff radius used by the neighbor list
    pub fn set_cutoff(&mut self, d0: f64) {
        self.neighbor_list_.set_cutoff(d0);
    }

    /// Provides the width of the buffer zone used by the neighbor list
    pub fn get_buffer_width(&self) -> f64 {
        self.neighbor_list_.buffer_width()
    }

    /// Modifies the width of the buffer zone used by the neighbor list
    pub fn set_buffer_width(&mut self, width: f64) {
        self.neighbor_list_.set_buffer_width(width);
    }

    /// Sets new rules that define which atoms and atom pairs can be neighbors.
    /// This method also triggers neighbor list update.
    pub fn set_rules(&mut self, nb_rules: Box<dyn NbListRules>) {
        self.neighbor_list_.set_rules(nb_rules);
        let stls_v = CoordinatesView {
            points: &self.coords_,
        };
        self.neighbor_list_.update_all_for_view(stls_v);
    }

    /// Calculates the periodic image of a vector relative to a reference point
    /*pub fn periodic_image(&self, reference: &Vec3, vector: &Vec3) -> Vec3 {
        let box_len = self.box_len();

        let mut image = vector.clone();
        let mut delta = reference.clone() - *vector;

        // Apply periodic boundary conditions along each dimension
        if delta.x > 0.5 * box_len {
            delta.x -= box_len;
            image.x += box_len;
        } else if delta.x < -0.5 * box_len {
            delta.x += box_len;
            image.x -= box_len;
        }

        if delta.y > 0.5 * box_len {
            delta.y -= box_len;
            image.y += box_len;
        } else if delta.y < -0.5 * box_len {
            delta.y += box_len;
            image.y -= box_len;
        }

        if delta.z > 0.5 * box_len {
            delta.z -= box_len;
            image.z += box_len;
        } else if delta.z < -0.5 * box_len {
            delta.z += box_len;
            image.z -= box_len;
        }

        image
    }*/
    /// Calculates the periodic image of a vector relative to a reference point
    pub fn get_periodic_image(&self, reference: &Vec3, vector: &Vec3) -> Vec3 {
        let box_len = self.get_box_len();

        let mut image = *vector;
        let delta = *reference - *vector;

        // Apply periodic boundary conditions along each dimension
        let half_box_len = 0.5 * box_len;

        // X-axis
        if delta.x.abs() > half_box_len {
            image.x += box_len * delta.x.signum();
        }

        // Y-axis
        if delta.y.abs() > half_box_len {
            image.y += box_len * delta.y.signum();
        }

        // Z-axis
        if delta.z.abs() > half_box_len {
            image.z += box_len * delta.z.signum();
        }
        image
    }
}

impl System for CartesianSystem {
    /// Returns the number of atoms of this system
    fn get_size(&self) -> usize {
        self.coords_.get_size()
    }

    /// Copies coordinates of a given atom from another system.
    ///
    /// This method also triggers neighbor list update at position `i`.
    /// # Arguments
    /// * `i` - index of an atom to be copied; it's the same index in both source and destination coordinates
    /// * `rhs` - the source to copy from
    fn copy_from(&mut self, i: usize, rhs: &Self) {
        self.coords_.copy_from(i, rhs.get_coordinates());
        let stls_v = CoordinatesView {
            points: &self.coords_,
        };
        self.neighbor_list_.update_for_view(stls_v, i);
    }
}

impl ResizableSystem for CartesianSystem {
    /// Change the number of atoms of this system
    fn set_size(&mut self, new_size: usize) {
        self.coords_.set_size(new_size);
    }

    /// Returns the maximum number of atoms system may have
    fn get_capacity(&self) -> usize {
        self.coords_.get_capacity()
    }
}
