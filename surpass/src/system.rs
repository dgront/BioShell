use std::collections::HashMap;
use std::fmt;
use std::fmt::Debug;
use num_traits::WrappingSub;
use bioshell_pdb::{load_pdb_file, Structure};
use bioshell_pdb::pdb_atom_filters::{IsCA, PdbAtomPredicate};
use bioshell_pdb::pdb_parsing_error::ParseError;

/// SURPASS-alpha system holds coordinates of all atoms
///
#[derive(Clone)]
pub struct SurpassAlphaSystem {
    chain_ids: Vec<u16>,
    cax: Vec<i32>,
    cay: Vec<i32>,
    caz: Vec<i32>,
    sgx: Vec<i32>,
    sgy: Vec<i32>,
    sgz: Vec<i32>,
    chains_map: HashMap<String, usize>,
    int_to_real: f64,
    int_to_real_2: f64
}

/// Converts real coordinate to its integer representation.
///
/// # Arguments
///  - `int_to_real_factor` - factor that converts integer coordinates into real ones in Angstroms
///  - `x` - the f64 coordinate value to be converted to i32
macro_rules! real_to_int {
    ($int_to_real_factor:expr, $x:expr) => {
        { ($x / $int_to_real_factor) as i32 }
    }
}

/// converts real coordinate to its integer representation.
///
/// # Arguments
///  - `int_to_real_factor` - factor that converts integer coordinates into real ones in Angstroms
///  - `ix` - the integer coordinate value to be converted to f64
macro_rules! int_to_real {
    ($int_to_real_factor:expr, $ix:expr) =>   {
        { $ix as f64 * $int_to_real_factor }
    }
}

impl SurpassAlphaSystem {

    fn new_empty(chain_lengths: &Vec<usize>, box_length: f64) -> SurpassAlphaSystem {
        let n_atoms = chain_lengths.iter().sum();
        let l = box_length / 2.0 / (i32::MAX as f64);
        let mut s = SurpassAlphaSystem{
            chain_ids: vec![0; n_atoms], cax: vec![0; n_atoms], cay: vec![0; n_atoms], caz: vec![0; n_atoms],
            sgx: vec![0; n_atoms], sgy: vec![0; n_atoms], sgz: vec![0; n_atoms],
            chains_map: Default::default(), int_to_real: l, int_to_real_2: l*l
        };
        let mut atoms_total = 0;
        let mut chain_idx = 0;
        for nc in chain_lengths {
            for i in 0..*nc {
                s.chain_ids[i+ atoms_total] = chain_idx;
            }
            chain_idx += 1;
            atoms_total += nc;
        }
        return s;
    }

    /// X coordinate of i-th CA atom
    pub fn cax(&self, i: usize) -> f64 { int_to_real!(self.int_to_real, self.cax[i]) }

    /// X coordinate of i-th CA atom
    pub fn cay(&self, i: usize) -> f64 { int_to_real!(self.int_to_real, self.cay[i]) }

    /// X coordinate of i-th CA atom
    pub fn caz(&self, i: usize) -> f64 { int_to_real!(self.int_to_real, self.caz[i]) }

    /// Index of a chain a given atom belongs to.
    ///
    /// # Arguments
    ///  -  `i` - atom index
    pub fn chain(&self, i: usize) -> u16 { self.chain_ids[i] }

    pub fn count_chains(&self) -> usize { self.chains_map.len() }

    /// Returns a string id denoting a given chain.
    ///
    /// # Arguments
    ///  -  `chain_idx` - integer index of a chain, starting from 0; to obtain the string id for
    ///     a particular atoms, use [`chain()`](chain()) method to find the integer index of that chain
    pub fn chain_id(&self, chain_idx: usize) -> Option<&String> {
        self.chains_map.iter()
            .find_map(|(key, &val)| if val == chain_idx { Some(key) } else { None })
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

        // pub fn create_random(nres: usize, nchains: usize) -> SurpassAlphaSystem {}

    pub fn from_pdb_structure(strctr: &Structure, box_length: f64) -> SurpassAlphaSystem {

        let mut chain_sizes: Vec<usize> = vec![];
        for chain_id in strctr.chain_ids().iter() {
            chain_sizes.push(strctr.chain_residue_ids(chain_id).len());
        }
        let mut surpass_model = SurpassAlphaSystem::new_empty(&chain_sizes, box_length);

        let is_ca = IsCA;
        let mut i_atom = 0;
        for ai in strctr.atoms().iter().filter(|&a| is_ca.check(a)) {
            if !surpass_model.chains_map.contains_key(&ai.chain_id) {
                surpass_model.chains_map.insert(ai.chain_id.clone(), surpass_model.chains_map.len());
            }
            surpass_model.chain_ids.push(*(surpass_model.chains_map.get(&ai.chain_id).unwrap()) as u16);
            surpass_model.cax[i_atom] = real_to_int!(surpass_model.int_to_real, ai.x);
            surpass_model.cay[i_atom] = real_to_int!(surpass_model.int_to_real, ai.y);
            surpass_model.caz[i_atom] = real_to_int!(surpass_model.int_to_real, ai.z);
            i_atom += 1;
        }
        assert_eq!(surpass_model.chains_map.len(), chain_sizes.len());

        return surpass_model;
    }

    pub fn from_pdb_file(fname: &str, box_length: f64) -> Result<SurpassAlphaSystem, ParseError> {
        
        let strctr = load_pdb_file(fname)?;
        return Ok(Self::from_pdb_structure(&strctr, box_length));
    }

    pub fn hinge_move(&mut self) {}
    pub fn tail_move(&mut self) {}
}


