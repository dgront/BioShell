use std::io::Write;
use std::ops::Range;
use rand::{Rng, thread_rng};
use bioshell_pdb::{load_pdb_file, Structure};
use bioshell_pdb::calc::Vec3;
use bioshell_pdb::nerf::restore_linear_chain;
use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
use bioshell_pdb::pdb_parsing_error::ParseError;
use bioshell_io::out_writer;

/// SURPASS-alpha system holds coordinates of all atoms
///
#[derive(Clone)]
pub struct SurpassAlphaSystem {
    pub cax: Vec<i32>,
    pub cay: Vec<i32>,
    pub caz: Vec<i32>,
    pub sgx: Vec<i32>,
    pub sgy: Vec<i32>,
    pub sgz: Vec<i32>,
    atoms_for_chain: Vec<Range<usize>>,
    chain_indexes: Vec<u16>,
    chain_id_by_index: Vec<String>,
    int_to_real: f64,
    int_to_real_2: f64,
    model_id: usize
}

impl SurpassAlphaSystem {

    pub fn new(chain_lengths: &[usize], box_length: f64) -> SurpassAlphaSystem {

        // ---------- Create an empty system
        let l = box_length / 2.0 / (i32::MAX as f64);
        let n_atoms = chain_lengths.iter().sum();
        let n_chains = chain_lengths.len();
        let mut s = SurpassAlphaSystem{
            chain_indexes: vec![0; n_atoms], cax: vec![0; n_atoms], cay: vec![0; n_atoms], caz: vec![0; n_atoms],
            sgx: vec![0; n_atoms], sgy: vec![0; n_atoms], sgz: vec![0; n_atoms],
            atoms_for_chain: vec![0..0; n_atoms], chain_id_by_index: vec![String::from("A");n_chains],
            int_to_real: l, int_to_real_2: l*l, model_id: 0
        };
        // ---------- Initialize coordinates
        let mut rnd = thread_rng();
        let r = vec![3.8; n_atoms];
        let mut planar: Vec<f64> = (0..n_atoms).map(|_| rnd.gen_range(90.0_f64.to_radians()..170.0_f64.to_radians())).collect();
        let mut dihedral: Vec<f64> = (0..n_atoms).map(|_| rnd.gen_range(-180.0_f64.to_radians()..190.0_f64.to_radians())).collect();
        let mut coords = vec![Vec3::default(); n_atoms];
        restore_linear_chain(&r[0..n_atoms], &planar[0..n_atoms], &dihedral[0..n_atoms], &mut coords[0..n_atoms]);
        for i in 0..n_atoms {
            s.vec3_to_ca(i,&coords[i]);
        }
        // ---------- Assign atoms to chains
        let mut atoms_total = 0;
        let codes: Vec<char> = ('A'..'Z').collect();
        for (ic,nc) in chain_lengths.iter().enumerate() {
            for i in 0..*nc { s.chain_indexes[i+ atoms_total] = ic as u16; }
            s.chain_id_by_index[ic] = codes[ic].to_string();
            s.atoms_for_chain[ic] = atoms_total..atoms_total+nc;
            atoms_total += nc;
        }
        return s;
    }

    #[inline(always)]
    pub fn int_to_real(&self, v: i32) -> f64 { self.int_to_real * v as f64 }

    #[inline(always)]
    pub fn real_to_int(&self, v: f64) -> i32 {
        (v / self.int_to_real).rem_euclid((i32::MAX as f64) + 1.0) as i32
    }

    /// Returns the number of atoms in this system (of all its chains)
    pub fn count_atoms(&self) -> usize { self.cax.len() }

    pub fn ca_to_vec3(&self, pos: usize) -> Vec3 {
        Vec3::new(self.int_to_real(self.cax[pos]),
                  self.int_to_real(self.cay[pos]),
                  self.int_to_real(self.caz[pos]))
    }

    pub fn vec3_to_ca(&mut self, pos: usize, v: &Vec3) {
        self.cax[pos] = self.real_to_int(v.x);
        self.cay[pos] = self.real_to_int(v.y);
        self.caz[pos] = self.real_to_int(v.z);
    }

    /// Index of a chain a given atom belongs to.
    ///
    /// # Arguments
    ///  -  `i` - atom index
    pub fn chain(&self, i: usize) -> u16 { self.chain_indexes[i] }

    /// Returns a range of atom indexes for a given chain of this system
    pub fn chain_atom(&self, ic: usize) -> &Range<usize> {
            &self.atoms_for_chain[ic]
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

        let mut r = self.cax[i].wrapping_sub(self.cax[j]);
        let mut r2 = r as f64 * r as f64;
        r = self.cay[i].wrapping_sub(self.cay[j]);
        r2 += r as f64 * r as f64;
        r = self.caz[i].wrapping_sub(self.caz[j]);
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
            surpass_model.cax[i_atom] = surpass_model.real_to_int(ai.pos.x);
            surpass_model.cay[i_atom] = surpass_model.real_to_int(ai.pos.y);
            surpass_model.caz[i_atom] = surpass_model.real_to_int(ai.pos.z);
            i_atom += 1;
        }
        assert_eq!(surpass_model.chain_id_by_index.len(), chain_sizes.len());

        return surpass_model;
    }

    pub fn from_pdb_file(fname: &str, box_length: f64) -> Result<SurpassAlphaSystem, ParseError> {
        
        let strctr = load_pdb_file(fname)?;
        return Ok(Self::from_pdb_structure(&strctr, box_length));
    }

    pub fn to_pdb_file(&self, fname: &str, if_append: bool) {

        let mut stream = out_writer(&fname, if_append);
        stream.write(format!("MODEL{:6}\n", self.model_id).as_bytes());
        for i in 0..self.count_atoms() {
            stream.write(format!("ATOM   {:4} {} GLY {}{:4}    {:8.3}{:8.3}{:8.3}  1.00  1.00           C\n",
                                 i, " CA ", &self.chain_id(self.chain(i) as usize).unwrap(), i,
                                 self.int_to_real(self.cax[i]),
                                 self.int_to_real(self.cay[i]),
                                 self.int_to_real(self.caz[i])).as_bytes()).expect("Error ocurred while writing PDB content!");

        }
        stream.write("ENDMDL\n".as_bytes());
    }
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

