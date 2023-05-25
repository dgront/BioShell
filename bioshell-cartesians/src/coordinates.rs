use std::ops::{Index, IndexMut, Range};

use bioshell_numerical::Vec3;
use bioshell_sim::{ResizableSystem, System};

/// Stateless immutable view of coordinates
pub struct CoordinatesView<'a> {
    pub points: &'a Coordinates,
}

/// Represents a system of Vec3 points in a simulation.
///
/// This system imposes simple cubic periodic boundary conditions on the modelled system
#[derive(Clone, Debug)]
pub struct Coordinates {
    box_len: f64,
    box_len_half: f64,
    current_size: usize,
    coords_vec: Vec<Vec3>,
    chains: Vec<Range<usize>>,
}

macro_rules! wrap_coordinate_to_box {
    ($val:expr, $L:expr, $coord:expr) => {
        $coord = $val;
        if $coord > $L {
            $coord = $coord - $L
        } else {
            if $coord < 0.0 {
                $coord = $L + $coord
            }
        }
    };
}

/// Calculates the shortest difference `c1 - c2`, taking periodic boundary condition into account
macro_rules! closest_image {
    ($c1:expr, $c2:expr, $L: expr,$L2: expr, $delta:expr) => {
        $delta = $c1 - $c2;
        if $delta > 0.0 {
            if $delta > $L2 {
                $delta -= $L
            }
        } else {
            if $delta < -$L2 {
                $delta += $L
            }
        }
    };
}

/// Finds the length of a simulation box to achieve assumed density
///
/// # Arguments
/// * `atom_radius` - atomic radius is same for all atoms of the system
/// * `n_atoms` - number of atoms in the system
/// * `density` - target density
pub fn box_width(atom_radius: f64, n_atoms: usize, density: f64) -> f64 {
    let v: f64 = 4.0 / 3.0 * std::f64::consts::PI * atom_radius.powi(3);
    (n_atoms as f64 * v / density).powf(1.0 / 3.0)
}

impl Coordinates {
    pub fn new(n: usize) -> Coordinates
    {
        let mut v = if n > 0
        {
            Vec::with_capacity(n)
        }
        else
        {
            Vec::new()
        };

        if n > 0
        {
            let zero = Vec3::from_float(0.0);
            v.resize(n, zero);
        }
        let chains: Vec<Range<usize>> = vec![0..n];
        let l: f64 = 100000.0;
        return Coordinates {
            current_size: 0,
            box_len: l,
            box_len_half: l / 2.0,
            coords_vec: v,
            chains,
        };
    }

    #[inline(always)]
    pub fn box_len(&self) -> f64 {
        self.box_len
    }

    #[inline(always)]
    pub fn set_box_len(&mut self, new_box_len: f64) {
        self.box_len = new_box_len;
        self.box_len_half = new_box_len / 2.0;
    }

    /// Returns the number of chains in this system
    pub fn count_chains(&self) -> usize {
        return self.chains.len();
    }

    /// Provides the (half-open) range of atoms that belong to a given chain.
    ///
    /// Per rust convention used in ``std::ops::Range`` struct, the returned
    /// ``start..end`` range contains all atoms indexed by ``start <= idx < end``
    pub fn chain_range(&self, idx: usize) -> &Range<usize> {
        &self.chains[idx]
    }

    /// Calculates the square of the (true) distance between two atoms
    ///
    /// # Arguments
    /// * `i` - index of the first atom
    /// * `j` - index of the second atom
    pub fn distance_square(&self, i: usize, j: usize) -> f64 {
        let mut d = self.coords_vec[i].x - self.coords_vec[j].x;
        let mut d2 = d * d;
        d = self.coords_vec[i].y - self.coords_vec[j].y;
        d2 += d * d;
        d = self.coords_vec[i].z - self.coords_vec[j].z;
        d2 += d * d;
        return d2;
    }

    /// Calculates the square of the shortest distance between two atoms.
    ///
    /// The returned distance is evaluate between atom `i` and its closest image of atom `j`
    ///
    /// # Arguments
    /// * `i` - index of the first atom
    /// * `j` - index of the second atom
    pub fn closest_distance_square(&self, i: usize, j: usize) -> f64 {
        let mut d: f64;
        closest_image!(self.coords_vec[i].x, self.coords_vec[j].x, self.box_len, self.box_len_half, d);
        let mut d2 = d * d;
        closest_image!(self.coords_vec[i].y, self.coords_vec[j].y, self.box_len, self.box_len_half, d);
        d2 += d * d;
        closest_image!(self.coords_vec[i].z, self.coords_vec[j].z, self.box_len, self.box_len_half, d);

        return d2 + d * d;
    }

    /// Calculates the square of the shortest distance between two atoms.
    ///
    /// The returned distance is evaluate between atom `i` and its closest image of a given vector
    ///
    /// # Arguments
    /// * `i` - index of the first atom
    /// * `v` - position of the second atom
    pub fn closest_distance_square_to_vec(&self, i: usize, v: &Vec3) -> f64 {
        let mut d: f64;
        closest_image!(self.coords_vec[i].x, v.x, self.box_len, self.box_len_half, d);
        let mut d2 = d * d;
        closest_image!(self.coords_vec[i].y, v.y, self.box_len, self.box_len_half, d);
        d2 += d * d;
        closest_image!(self.coords_vec[i].z, v.z, self.box_len, self.box_len_half, d);

        return d2 + d * d;
    }

    /// Creates a Vec3 that holds the image of `the_atom` that is the closest to `ref_atom`
    ///
    /// # Arguments
    /// * `ref_atom` - the reference atom
    /// * `the_atom` - index of the atom to be cloned
    pub fn clone_closest_image(&self, ref_atom: usize, the_atom: usize) -> Vec3 {

        let mut out = self.coords_vec[ref_atom].clone();
        let mut d: f64;
        closest_image!(self.coords_vec[the_atom].x, self.coords_vec[ref_atom].x, self.box_len, self.box_len_half, d);
        out.x += d;
        closest_image!(self.coords_vec[the_atom].y, self.coords_vec[ref_atom].y, self.box_len, self.box_len_half, d);
        out.y += d;
        closest_image!(self.coords_vec[the_atom].z, self.coords_vec[ref_atom].z, self.box_len, self.box_len_half, d);
        out.z += d;

        return out;
    }

    /// Calculates the difference in ``x`` coordinate between the i-th atom and a given ``x`` value
    /// This function obeys periodic boundary conditions and returns the distance to the closest
    /// image of the  position ``i``
    pub fn delta_x(&self, i: usize, x: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.coords_vec[i].x, x, self.box_len, self.box_len_half, d);
        d
    }

    /// Calculates the difference in ``y`` coordinate between the i-th atom and a given ``y`` value
    /// This function obeys periodic boundary conditions and returns the distance to the closest
    /// image of the  position ``i``
    pub fn delta_y(&self, i: usize, y: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.coords_vec[i].y, y, self.box_len, self.box_len_half, d);
        d
    }

    /// Calculates the difference in ``z`` coordinate between the i-th atom and a given ``z`` value
    /// This function obeys periodic boundary conditions and returns the distance to the closest
    /// image of the  position ``i``
    pub fn delta_z(&self, i: usize, z: f64) -> f64 {
        let mut d: f64;
        closest_image!(self.coords_vec[i].z, z, self.box_len, self.box_len_half, d);
        d
    }

    /// Provides the chain index type of the `i`-th atom of this system
    pub fn chain_id(&self, i: usize) -> u16 {
        self.coords_vec[i].chain_id
    }

    /// Provides the residue type of the `i`-th **atom** of this system
    pub fn res_type(&self, i: usize) -> u8 {
        self.coords_vec[i].res_type
    }

    /// Provides the type of the `i`-th atom of this system
    pub fn atom_type(&self, i: usize) -> u8 {
        self.coords_vec[i].atom_type
    }

    /// Assign each atom of this system to a chain
    ///
    /// This method assigns a new chain index for each atom based of a vector of atom ranges.
    /// The number of chains of this system will we equal to `chain_ranges.size()`,
    /// the `k`-th chain will include atoms `$r0 \le i \t r1$` where `$r0$` and `$r1$` are defined
    /// by the `k`-th range: from `chain_ranges[k].0` to `chain_ranges[k].1`
    ///
    /// Assigning all the chains for this [`Coordinates`] help avoiding incorrect assignments
    pub fn set_chains(&mut self, chain_ranges: &Vec<(usize, usize)>) {
        assert_eq!(
            chain_ranges.last().unwrap().1,
            self.current_size,
            "chain ranges do not mach the size of the coordinates"
        );
        let mut id: u16 = 0;
        for r in chain_ranges
        {
            for i in r.0..r.1 {
                self.coords_vec[i].chain_id = id;
            }
            id += 1;
        }
    }

    /// Assigns the residue type for the **atom** `i` to `t`
    pub fn set_res_type(&mut self, i: usize, t: u8) {
        self.coords_vec[i].res_type = t;
    }

    /// Assigns the type of the atom `i` to `t`
    pub fn set_atom_type(&mut self, i: usize, t: u8) {
        self.coords_vec[i].atom_type = t;
    }

    /// 'x' coordinate of `i`-th atom of this system
    pub fn x(&self, i: usize) -> f64 {
        self.coords_vec[i].x
    }

    /// 'y' coordinate of `i`-th atom of this system
    pub fn y(&self, i: usize) -> f64 {
        self.coords_vec[i].y
    }

    /// 'z' coordinate of `i`-th atom of this system
    pub fn z(&self, i: usize) -> f64 {
        self.coords_vec[i].z
    }

    pub fn set_x(&mut self, i: usize, x: f64) {
        wrap_coordinate_to_box!(x, self.box_len, self.coords_vec[i].x);
    }

    pub fn set_y(&mut self, i: usize, y: f64) {
        wrap_coordinate_to_box!(y, self.box_len, self.coords_vec[i].y);
    }

    pub fn set_z(&mut self, i: usize, z: f64) {
        wrap_coordinate_to_box!(z, self.box_len, self.coords_vec[i].z);
    }

    pub fn set(&mut self, i: usize, x: f64, y: f64, z: f64) {
        wrap_coordinate_to_box!(x, self.box_len, self.coords_vec[i].x);
        wrap_coordinate_to_box!(y, self.box_len, self.coords_vec[i].y);
        wrap_coordinate_to_box!(z, self.box_len, self.coords_vec[i].z);
    }

    pub fn set_vec(&mut self, i: usize, vec: Vec3) {
        let x = vec.x;
        let y = vec.y;
        let z = vec.z;
        self.set(i, x, y, z);
    }

    pub fn add(&mut self, i: usize, x: f64, y: f64, z: f64) {
        wrap_coordinate_to_box!(self.coords_vec[i].x + x, self.box_len, self.coords_vec[i].x);
        wrap_coordinate_to_box!(self.coords_vec[i].y + y, self.box_len, self.coords_vec[i].y);
        wrap_coordinate_to_box!(self.coords_vec[i].z + z, self.box_len, self.coords_vec[i].z);
    }

    /// Copy coordinates of i-th atom from a given vector.
    ///
    /// This method (unlike the [`set()`](set) method) does not apply PBC. To the contrary,
    /// it assumes the two systems: `self` and `rhs` have exactly the same simulation box geometry
    pub fn copy_from_vec(&mut self, i: usize, rhs: &Vec3) {
        self.coords_vec[i].x = rhs.x;
        self.coords_vec[i].y = rhs.y;
        self.coords_vec[i].z = rhs.z;
    }
}

impl Index<usize> for Coordinates {
    type Output = Vec3;
    fn index(&self, i: usize) -> &Vec3 {
        &self.coords_vec[i]
    }
}

impl IndexMut<usize> for Coordinates {
    fn index_mut(&mut self, i: usize) -> &mut Vec3 {
        &mut self.coords_vec[i]
    }
}

impl System for Coordinates {
    /// Returns the current number of atoms of this system
    fn size(&self) -> usize {
        return self.current_size;
    }

    /// Copy coordinates of i-th atom from a given `rhs` coordinates.
    ///
    /// This method (unlike the [`set()`](set) method) does not apply PBC. To the contrary,
    /// it assumes the two systems: `self` and `rhs` have exactly the same simulation box geometry
    fn copy_from(&mut self, i: usize, rhs: &Self) {
        self.coords_vec[i].x = rhs.coords_vec[i].x;
        self.coords_vec[i].y = rhs.coords_vec[i].y;
        self.coords_vec[i].z = rhs.coords_vec[i].z;
    }
}

impl ResizableSystem for Coordinates {
    /// Changes the number of atoms of this system
    fn set_size(&mut self, new_size: usize) {
        self.current_size = new_size;
    }

    /// Returns the maximum number of atoms of this system
    fn capacity(&self) -> usize {
        return self.coords_vec.len();
    }
}
