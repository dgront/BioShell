use std::io::Write;
use std::ops::Range;
use rand::{Rng};
use bioshell_pdb::{load_pdb_file, Structure};
use bioshell_pdb::calc::Vec3;
use bioshell_pdb::nerf::restore_linear_chain;
use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
use bioshell_pdb::PDBError;
use bioshell_io::out_writer;


/// SURPASS-alpha system holds coordinates of all atoms
///
#[derive(Clone)]
pub struct SurpassAlphaSystem {
    pub bbx: Vec<i32>,
    pub bby: Vec<i32>,
    pub bbz: Vec<i32>,
    pub sgx: Vec<i32>,
    pub sgy: Vec<i32>,
    pub sgz: Vec<i32>,
    residues_for_chain: Vec<Range<usize>>,
    chain_indexes: Vec<u16>,
    chain_id_by_index: Vec<String>,
    int_to_real: f64,
    int_to_real_2: f64,
    box_length: f64,
    model_id: usize
}

macro_rules! closest_coordinate {
    ($self: expr, $cvec: expr, $pos: expr, $ref_pos: expr) => {
        ($cvec[$pos].wrapping_sub($cvec[$ref_pos]) as f64 + $cvec[$ref_pos] as f64) * $self.int_to_real
    }
}
impl SurpassAlphaSystem {

    pub fn new(chain_lengths: &[usize], box_length: f64) -> SurpassAlphaSystem {

        // ---------- Create an empty system
        let l = box_length / 2.0 / (i32::MAX as f64);
        let n_atoms = chain_lengths.iter().sum();
        let n_chains = chain_lengths.len();
        let mut s = SurpassAlphaSystem{
            chain_indexes: vec![0; n_atoms], bbx: vec![0; n_atoms], bby: vec![0; n_atoms], bbz: vec![0; n_atoms],
            sgx: vec![0; n_atoms], sgy: vec![0; n_atoms], sgz: vec![0; n_atoms],
            residues_for_chain: vec![0..0; n_atoms], chain_id_by_index: vec![String::from("A"); n_chains],
            int_to_real: l, int_to_real_2: l*l, box_length, model_id: 0
        };
        // ---------- Assign atoms to chains
        let mut atoms_total = 0;
        let codes: Vec<char> = ('A'..'Z').collect();
        for (ic,nc) in chain_lengths.iter().enumerate() {
            for i in 0..*nc { s.chain_indexes[i+ atoms_total] = ic as u16; }
            s.chain_id_by_index[ic] = codes[ic].to_string();
            s.residues_for_chain[ic] = atoms_total..atoms_total+nc;
            atoms_total += nc;
        }

        return s;
    }

    #[inline(always)]
    pub fn int_to_real(&self, v: i32) -> f64 { self.int_to_real * v as f64 }

    #[inline(always)]
    pub fn real_to_int(&self, v: f64) -> i32 {
        (v / self.int_to_real).rem_euclid((u32::MAX as f64) + 1.0)  as u32 as i32
    }

    /// The length of a simulation box
    pub fn box_length(&self) -> f64 { self.box_length }

    /// Adjusts all bonds for the given length.
    ///
    /// Every bond will get shortened or extended along ist direction to match the given new length.
    /// This procedure doesn't alter any planar or dihedral angle. Every chain is treated separately,
    /// which means that every first atom of each chain will remain at its original position.
    /// A whole chain may be longer or shorter, depending on the circumstances.
    ///
    /// The intended use of this method is to correct bond lengths which after numerous rotational
    /// Monte Carlo moves may slightly diverge from their assumed length.
    pub fn adjust_bond_length(&mut self, new_length: f64) -> f64 {
        let mut mav_violation = 0.0;
        for i_chain in 0..self.count_chains() {
            let from_to = self.residues_for_chain[i_chain].clone();
            let mut old_previous = self.atom_to_vec3(from_to.start);
            let mut updated_previous = self.atom_to_vec3(from_to.start);
            let mut current = Vec3::default();
            for i_res in from_to.start+1..from_to.end {
                self.set_ca_to_nearest_vec3(i_res, i_res-1, &mut current);
                current -= &old_previous;
                mav_violation = (current.length()-new_length).abs().max(mav_violation);
                current.normalize();
                current *= new_length;
                current += &updated_previous;
                self.set_ca_to_nearest_vec3(i_res, i_res-1, &mut old_previous);
                self.vec3_to_ca(i_res, &current);
                updated_previous.set(&current);
            }
        }

        return mav_violation;
    }

    /// Returns the number of atoms in this system (of all its chains)
    pub fn count_residues(&self) -> usize { self.bbx.len() }

    /// Returns real coordinates of a C-alpha atom in a newly created Vec3 struct.
    ///
    /// Coordinates of a `pos` alpha carbon are converted from their integer representation
    /// to `f64` values and returned as a [`Vec3`](bioshell_pdb::calc::vec3) struct.
    pub fn atom_to_vec3(&self, pos: usize) -> Vec3 {
        Vec3::new(self.int_to_real(self.bbx[pos]),
                  self.int_to_real(self.bby[pos]),
                  self.int_to_real(self.bbz[pos]))
    }
    /// Copies real coordinates of a C-alpha atom into a given Vec3 struct.
    ///
    /// Coordinates of a `pos` alpha carbon are converted from their integer representation
    /// to `f64` values and stored in a given [`Vec3`](bioshell_pdb::calc::vec3) struct.
    pub fn set_atom_to_vec3(&self, pos: usize, v: &mut Vec3) {
        v.set3(self.int_to_real(self.bbx[pos]),
               self.int_to_real(self.bby[pos]),
               self.int_to_real(self.bbz[pos]));
    }

    /// Returns real coordinates of an image of a C-alpha atom that is the closest to a given atom.
    ///
    /// This method finds the periodic image of a `pos` C-alpha atom that is the closes to the
    /// `ref_pos` C-alpha atom  and returns its coordinates as a newly created  [`Vec3`](bioshell_pdb::calc::vec3) struct.
    pub fn atom_to_nearest_vec3(&self, pos: usize, ref_pos: usize) -> Vec3 {

        Vec3::new(closest_coordinate!(self, self.bbx, pos, ref_pos),
                  closest_coordinate!(self, self.bby, pos, ref_pos),
                  closest_coordinate!(self, self.bbz, pos, ref_pos))
    }

    /// Returns real coordinates of an image of a C-alpha atom that is the closest to a given atom.
    ///
    /// This method finds the periodic image of a `pos` C-alpha atom that is the closes to the
    /// `ref_pos` C-alpha atom  and stores its coordinates in a given  [`Vec3`](bioshell_pdb::calc::vec3) struct.
    pub fn set_ca_to_nearest_vec3(&self, pos: usize, ref_pos: usize, v: &mut Vec3) {

        v.set3(closest_coordinate!(self, self.bbx, pos, ref_pos),
               closest_coordinate!(self, self.bby, pos, ref_pos),
               closest_coordinate!(self, self.bbz, pos, ref_pos));
    }

    pub fn vec3_to_ca(&mut self, pos: usize, v: &Vec3) {
        self.bbx[pos] = self.real_to_int(v.x);
        self.bby[pos] = self.real_to_int(v.y);
        self.bbz[pos] = self.real_to_int(v.z);
    }

    /// Index of a chain a given atom belongs to.
    ///
    /// # Arguments
    ///  -  `i` - atom index
    pub fn chain(&self, i: usize) -> u16 { self.chain_indexes[i] }

    /// Returns a range of atom indexes for a given chain of this system
    pub fn chain_residues(&self, ic: usize) -> &Range<usize> {
            &self.residues_for_chain[ic]
    }

    /// Returns the number of chains in this system
    pub fn count_chains(&self) -> usize { self.chain_id_by_index.len() }

    /// Returns a string id denoting a given chain.
    ///
    /// # Arguments
    ///  -  `chain_idx` - integer index of a chain, starting from 0; to obtain the string id for
    ///     a particular atoms, use [`chain()`](chain()) method to find the integer index of that chain
    pub fn chain_id(&self, chain_idx: usize) -> Result<&String, &str> {

        if chain_idx < self.chain_id_by_index.len() { return Ok(&self.chain_id_by_index[chain_idx])}
        else {return Err("Incorrect chain index") }
    }

    /// Set a new chain ID string for one of the chains of this system.
    ///
    /// ```
    /// # use surpass::SurpassAlphaSystem;
    /// # let mut system = SurpassAlphaSystem::new(&[5, 5], 100.0);
    /// // by default the first chain is named as "A"
    /// assert_eq!(system.chain_id(0).unwrap(), &String::from("A"));
    /// system.set_chain_id(0, "XYZ");
    /// assert_eq!(system.chain_id(0).unwrap(), &String::from("XYZ"));
    /// system.set_chain_id(0, "ZYX");
    /// assert_eq!(system.chain_id(0).unwrap(), &String::from("ZYX"));
    /// ```
    pub fn set_chain_id(&mut self, chain_idx: usize, chain_id: &str) {
        self.chain_id_by_index[chain_idx] = String::from(chain_id);
    }

    pub fn distance_squared(&self, i: usize, j: usize) -> f64 {

        let mut r = self.bbx[i].wrapping_sub(self.bbx[j]);
        let mut r2 = r as f64 * r as f64;
        r = self.bby[i].wrapping_sub(self.bby[j]);
        r2 += r as f64 * r as f64;
        r = self.bbz[i].wrapping_sub(self.bbz[j]);
        r2 += r as f64 * r as f64;

        return r2 * self.int_to_real_2;
    }

    pub fn distance(&self, i: usize, j: usize) -> f64 { self.distance_squared(i, j).sqrt() }

    pub fn from_pdb_structure(strctr: &Structure, box_length: f64) -> SurpassAlphaSystem {

        let mut chain_sizes: Vec<usize> = vec![];
        for chain_id in strctr.chain_ids().iter() {
            chain_sizes.push(strctr.chain_residue_ids(chain_id).len());
        }
        let mut surpass_model = SurpassAlphaSystem::new(&chain_sizes, box_length);
        for (i, cid) in strctr.chain_ids().iter().enumerate() {
            surpass_model.set_chain_id(i, cid);
        }

        let is_ca = IsCA;
        let mut i_atom = 0;
        for ai in strctr.atoms().iter().filter(|&a| is_ca.check(a)) {
            surpass_model.bbx[i_atom] = surpass_model.real_to_int(ai.pos.x);
            surpass_model.bby[i_atom] = surpass_model.real_to_int(ai.pos.y);
            surpass_model.bbz[i_atom] = surpass_model.real_to_int(ai.pos.z);
            i_atom += 1;
        }
        assert_eq!(surpass_model.chain_id_by_index.len(), chain_sizes.len());

        return surpass_model;
    }

    pub fn make_random<R: Rng>(chain_lengths: &[usize], box_length: f64, rnd_gen: &mut R) -> SurpassAlphaSystem {
        let mut s = SurpassAlphaSystem::new(chain_lengths, box_length);
        let n_atoms = s.count_residues();
        // ---------- Initialize coordinates
        let r = vec![3.8; n_atoms];
        let planar: Vec<f64> = (0..n_atoms).map(|_| rnd_gen.gen_range(150.0_f64.to_radians()..170.0_f64.to_radians())).collect();
        let dihedral: Vec<f64> = (0..n_atoms).map(|_| rnd_gen.gen_range(-180.0_f64.to_radians()..190.0_f64.to_radians())).collect();
        let mut coords = vec![Vec3::default(); n_atoms];
        restore_linear_chain(&r[0..n_atoms], &planar[0..n_atoms], &dihedral[0..n_atoms], &mut coords[0..n_atoms]);
        for i in 0..n_atoms {
            s.vec3_to_ca(i,&coords[i]);
        }

        return s;
    }

    pub fn from_pdb_file(fname: &str, box_length: f64) -> Result<SurpassAlphaSystem, PDBError> {
        
        let strctr = load_pdb_file(fname)?;
        return Ok(Self::from_pdb_structure(&strctr, box_length));
    }

    pub fn to_pdb_file(&self, fname: &str, if_append: bool) {

        let mut stream = out_writer(&fname, if_append);
        stream.write(format!("MODEL{:6}\n", self.model_id).as_bytes());
        for i_chain in 0..self.count_chains() {
            let begin_atom_idx = self.residues_for_chain[i_chain].start;
            let end_atom_idx = self.residues_for_chain[i_chain].end;
            let chain_code = &self.chain_id_by_index[i_chain];
            let mut the_atom: Vec3 = self.atom_to_vec3(begin_atom_idx);
            let mut i_resid = 1;
            stream.write(format!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00  1.00           C\n",
                                 begin_atom_idx, " CA ", chain_code, i_resid,
                                 the_atom.x,
                                 the_atom.y,
                                 the_atom.z).as_bytes()).expect("Error occurred while writing PDB content!");
            for i_atom in begin_atom_idx+1..end_atom_idx {
                i_resid += 1;
                self.set_ca_to_nearest_vec3(i_atom, begin_atom_idx, &mut the_atom);
                stream.write(format!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  2.00  2.00           C\n",
                                     i_atom, " CA ", chain_code, i_resid,
                                     the_atom.x,
                                     the_atom.y,
                                     the_atom.z).as_bytes()).expect("Error occurred while writing PDB content!");
            }

        }
        stream.write("ENDMDL\n".as_bytes());
    }
}

pub fn calculate_cm(system: &SurpassAlphaSystem, i_chain: usize) -> Vec3 {

    let begin_atom_idx = system.residues_for_chain[i_chain].start;
    let end_atom_idx = system.residues_for_chain[i_chain].end;
    let mut cm_vec: Vec3 = system.atom_to_vec3(begin_atom_idx);
    let mut tmp_atom: Vec3 = system.atom_to_vec3(begin_atom_idx);
    for i_atom in begin_atom_idx+1..end_atom_idx {
        system.set_ca_to_nearest_vec3(i_atom, begin_atom_idx, &mut tmp_atom);
        cm_vec += &tmp_atom;
    }
    cm_vec /= system.chain_residues(i_chain).len() as f64;

    return cm_vec;
}

/// Creates a system that contains a single chain in an extended conformation
///
/// The newly created system will contain a single chain of `n_res` residues, placed in the middle of a cubic periodic box
/// of length `box_length`. All planar and dihedral angles are set to 120 and 180 degrees, respectively
pub fn extended_chain(n_res: usize, box_length: f64) -> SurpassAlphaSystem {

    let mut model = SurpassAlphaSystem::new(&[n_res], box_length);
    // ---------- Initialize internal coordinates
    let r= vec![3.8; n_res];
    let planar = vec![120.0_f64.to_radians(); n_res];
    let dihedral = vec![180.0_f64.to_radians(); n_res];
    let mut coords = vec![Vec3::default(); n_res];
    restore_linear_chain(&r, &planar, &dihedral, &mut coords);
    for i in 0..n_res {
        model.vec3_to_ca(i, &coords[i]);
    }

    return model;
}

