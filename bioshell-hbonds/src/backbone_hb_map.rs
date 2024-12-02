use std::collections::HashMap;
use bioshell_pdb::{PdbAtom, ResidueId, Structure};
use bioshell_pdb::calc::{planar_angle3, Vec3};
use crate::{MAX_AD_DISTANCE, MAX_AH_DISTANCE, MIN_AHD_ANGLE, MIN_PAH_ANGLE, N_H_BOND_LENGTH};

const DSSP_CONST: f64  = 0.42 * 0.2 * 332.0;

/// Provides an integer index for each residue in a structure.
///
/// The struct converts a residue id to an integer index and vice versa.
pub struct ResidueIndexer {
    residue_index: HashMap<ResidueId, usize>,
}

impl ResidueIndexer {

    pub fn from_structure(strctr: &Structure) -> ResidueIndexer {
        let mut residue_index: HashMap<ResidueId, usize> = strctr.residue_ids().iter()
            .filter(|id| strctr.residue_type(id).unwrap().chem_compound_type.is_peptide_linking())
            .enumerate()
            .map(|(i, id)| (id.clone(), i))
            .collect();

        ResidueIndexer { residue_index }
    }

    pub fn index(&self, residue_id: &ResidueId) -> Option<usize> {
        self.residue_index.get(residue_id).copied()
    }

    pub fn residue_id(&self, index: usize) -> &ResidueId {
        self.residue_index.iter().find(|(_, i)| **i == index).unwrap().0
    }

    pub fn len(&self) -> usize { self.residue_index.len() }
}



pub struct BackboneHBond<'a> {
    the_structure: &'a Structure,
    donated_h: Vec3,
    donor_n: usize,
    acceptor_o: usize,
    acceptor_c: usize,
}

impl<'a> BackboneHBond<'a> {

    pub fn new(strctr: &'a Structure, donated_h: Vec3,
                donor_n: usize, acceptor_o: usize, acceptor_c: usize) -> BackboneHBond<'a> {

        BackboneHBond{ the_structure: strctr, donated_h, donor_n, acceptor_o, acceptor_c, }
    }

    #[allow(non_snake_case)]
    pub fn distance_AH(&self) -> f64 { self.the_structure.atoms()[self.acceptor_o].pos.distance_to(&self.donated_h) }

    #[allow(non_snake_case)]
    pub fn distance_DA(&self) -> f64 {
        let n_atom = &self.the_structure.atoms()[self.donor_n];
        self.the_structure.atoms()[self.acceptor_o].pos.distance_to(&n_atom.pos)
    }

    #[allow(non_snake_case)]
    pub fn angle_AHD(&self) -> f64 {
        let n_atom = &self.the_structure.atoms()[self.donor_n];
        let o_atom = &self.the_structure.atoms()[self.acceptor_o];
        planar_angle3(&n_atom.pos, &self.donated_h, &o_atom.pos)
    }

    #[allow(non_snake_case)]
    pub fn angle_PAH(&self) -> f64 {
        let c_atom = &self.the_structure.atoms()[self.acceptor_c];
        let o_atom = &self.the_structure.atoms()[self.acceptor_o];
        planar_angle3(&c_atom.pos, &o_atom.pos, &self.donated_h)
    }

    /// Calculates the DSSP energy for a bond.
    ///
    /// DSSP energy was defined by Kabsch and Sander using the following formula:
    ///
    /// ```math
    /// E = 0.084 { 1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN } * 332 kcal/mol
    /// ```
    ///
    /// # References
    /// Kabsch W, Sander C (1983). "Dictionary of protein secondary structure: pattern recognition of
    /// hydrogen-bonded and geometrical features". Biopolymers 22 (12): 2577–637.
    pub fn dssp_energy(&self) -> f64 {
        let donor_atom = &self.the_structure.atoms()[self.donor_n].pos;
        let hydrogen = &self.donated_h;
        let acceptor_atom = &self.the_structure.atoms()[self.acceptor_o].pos;
        let behind_acceptor_atom = &self.the_structure.atoms()[self.acceptor_c].pos;
        let mut e = 1.0 / donor_atom.distance_to(acceptor_atom); // N-O
        e += 1.0 / hydrogen.distance_to(behind_acceptor_atom); // C-H
        e -= 1.0 / hydrogen.distance_to(acceptor_atom); // O-H
        e -= 1.0 / donor_atom.distance_to(behind_acceptor_atom); // N-C

        return e * DSSP_CONST;
    }
}

pub struct BackboneHBondMap<'a> {
    the_structure: &'a Structure,
    hydrogens: Vec<PdbAtom>,
    indexer: ResidueIndexer,
    h_bonds: HashMap<(usize,usize), BackboneHBond<'a>>,
}

impl<'a> BackboneHBondMap<'a> {

    pub fn new(strctr: &'a Structure) -> Self {
        let mut hbonds = BackboneHBondMap {
            the_structure: strctr,
            hydrogens: vec![],
            indexer: ResidueIndexer::from_structure(strctr),
            h_bonds: Default::default()
        };

        hbonds.find_donors_acceptors();

        return hbonds;
    }

    pub fn h_bonds(&self)-> impl Iterator<Item = (&(usize, usize), &BackboneHBond<'a>)> {
        self.h_bonds.iter()
    }
    pub fn get_h_bond(&self, donor_residue: &ResidueId, acceptor_residue: &ResidueId) -> Option<&BackboneHBond<'a>> {

        let d = self.indexer.index(donor_residue);
        let a = self.indexer.index(acceptor_residue);
        if let (Some(d), Some(a)) = (d, a) {
            self.h_bonds.get(&(d, a));
        } else { return None }

        return return None;
    }

    fn find_donors_acceptors(&mut self) {

        // ---------- Extract backbone atoms ------------
        let mut bb: Vec<BackboneResidue> = vec![BackboneResidue::new(); self.indexer.len()];
        for (idx, atom) in self.the_structure.atoms().iter().enumerate() {
            let res_id = ResidueId::try_from(atom).unwrap();
            if let Some(res_idx) = self.indexer.index(&res_id) {
                match atom.name.as_str() {
                    " N  " => { bb[res_idx].n = Some(idx); }
                    " CA " => { bb[res_idx].ca = Some(idx); }
                    " C  " => { bb[res_idx].c = Some(idx); }
                    " O  " => { bb[res_idx].o = Some(idx); }
                    _ => {}
                }
            }
        }

        // ---------- Calculate positions of hydrogen atoms ------------
        let mut hydrogens: Vec<Option<Vec3>> = vec![None; self.indexer.len()];
        for (idx, res) in bb.iter().enumerate() {
            if !res.is_complete() { continue; }
            if idx == 0 { continue; }
            if let Some(prev_c) = bb[idx-1].c {
                let n = res.n.unwrap();
                let ca = res.ca.unwrap();
                let n_pos = self.the_structure.atoms()[n].pos;
                let ca_pos = self.the_structure.atoms()[ca].pos;
                let c_pos = self.the_structure.atoms()[prev_c].pos;
                let h = peptide_hydrogen(&c_pos, &n_pos, &ca_pos);
                hydrogens[idx] = Some(h);
            }
        }

        // ---------- Detect hydrogen bonds ------------
        let atoms = self.the_structure.atoms(); // --- alias to avoid repeated calls
        for (d_idx, d_res) in bb.iter().enumerate() {
            if !d_res.is_complete() { continue; }
            if hydrogens[d_idx].is_none() { continue; }
            let d_n = d_res.n.unwrap();
            for (a_idx, a_res) in bb.iter().enumerate() {
                if !a_res.is_complete() { continue; }
                let a_o = a_res.o.unwrap();
                let a_c = a_res.c.unwrap();
                // --- create the H-bond object
                let hb = BackboneHBond::new(self.the_structure, hydrogens[d_idx].unwrap().clone(), d_n, a_o, a_c);
                // --- check distances and angles
                if hb.distance_AH() > MAX_AH_DISTANCE { continue; }
                if hb.distance_DA() > MAX_AD_DISTANCE { continue; }
                if hb.angle_AHD().to_degrees() < MIN_AHD_ANGLE { continue; }
                if hb.angle_PAH().to_degrees() < MIN_PAH_ANGLE { continue; }
                // --- it's OK, insert that into the map
                self.h_bonds.insert((d_idx, a_idx), hb);
            }
        }
    }
}

pub fn peptide_hydrogen(prev_c: &Vec3, the_n: &Vec3, the_ca: &Vec3) -> Vec3 {

    let mut h = Vec3::new(the_n.x, the_n.y, the_n.z);
    h -= prev_c;
    h.normalize();

    let mut tmp2 = Vec3::new(the_n.x, the_n.y, the_n.z);
    tmp2 -= the_ca;
    tmp2.normalize();

    h += &tmp2;
    h.normalize();
    h *= N_H_BOND_LENGTH;
    h += &the_n;

    return h;
}

#[derive(Clone)]
struct BackboneResidue {
    n: Option<usize>,
    ca: Option<usize>,
    c: Option<usize>,
    o: Option<usize>,
}

impl BackboneResidue {
    pub fn new() -> BackboneResidue { BackboneResidue { n: None, ca: None, c: None, o: None} }

    pub fn is_complete(&self) -> bool {
        self.n.is_some() && self.ca.is_some() && self.c.is_some() && self.o.is_some()
    }
}