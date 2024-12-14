use std::collections::HashMap;
use bioshell_pdb::{PdbAtom, ResidueId, Structure};
use bioshell_pdb::calc::{planar_angle3, Vec3};
use crate::{MAX_AD_DISTANCE, MAX_AH_DISTANCE, MIN_AHD_ANGLE, MIN_PAH_ANGLE, N_H_BOND_LENGTH};

const DSSP_CONST: f64  = 0.42 * 0.2 * 332.0;

/// Provides an integer index for each residue in a structure.
///
/// The struct converts a residue id to an integer index and vice versa.
pub struct ResidueIndexer {
    residue_ids: Vec<ResidueId>,
    residue_index: HashMap<ResidueId, usize>,
}

impl ResidueIndexer {

    pub fn from_structure(strctr: &Structure) -> ResidueIndexer {

        let mut residue_ids: Vec<ResidueId> = vec![];
        let mut residue_index: HashMap<ResidueId, usize> = Default::default();

        let mut idx = 0;
        for residue_id in strctr.residue_ids() {
            if !strctr.residue_type(residue_id).unwrap().chem_compound_type.is_peptide_linking() { continue; }
            residue_index.insert(residue_id.clone(), idx);
            residue_ids.push(residue_id.clone());
            idx += 1;
        }
        ResidueIndexer { residue_ids, residue_index }
    }

    /// Returns the index of the given residue.
    ///
    /// The opposite operation to this is [`ResidueIndexer::residue_id()`].
    pub fn index(&self, residue_id: &ResidueId) -> Option<usize> { self.residue_index.get(residue_id).copied() }

    /// Iterates over all the residues in the structure indexed by this [`ResidueIndexer`].
    pub fn residue_ids(&self) -> impl Iterator<Item=&ResidueId> { self.residue_ids.iter() }

    /// Returns the [`ResidueId`] for the given residue index.
    ///
    /// The opposite operation to this is [`ResidueIndexer::index()`].
    pub fn residue_id(&self, index: usize) -> &ResidueId { &self.residue_ids[index] }

    /// Provides the number of residues in the structure indexed by this object.
    ///
    /// Note, that not all the residues of the structure may be indexed. Some, like ions and water molecules
    /// are omitted.
    pub fn n_residues(&self) -> usize { self.residue_index.len() }
}



/// A single hydrogen bond between a backbone nitrogen and an oxygen atoms.
///
/// [`BackboneHBond`] provides a detailed description of a protein backbone hydrogen bond:
/// ```text
///     O=C             N-H
///        \           /
///         N-H ... O=C
///        /           \
///     R-C             C-R
/// ```
/// marked as `N-H ... O=C`, between an amide group `N-H` donating its  hydrogen `H` to a carbonyl group `O=C` accepting the hydrogen.
/// The `N` atom is referred to as the `donor` (`'D'`) and the `O` atom is referred to as the `acceptor` (`'A'`).
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

    /// Calculates the distance between the hydrogen atom, and it's acceptor (a carbonyl oxygen in the case of a backbone hydrogen bond).
    #[allow(non_snake_case)]
    pub fn distance_AH(&self) -> f64 { self.the_structure.atoms()[self.acceptor_o].pos.distance_to(&self.donated_h) }

    /// Calculates the distance between the donor and acceptor atoms
    #[allow(non_snake_case)]
    pub fn distance_DA(&self) -> f64 {
        let n_atom = &self.the_structure.atoms()[self.donor_n];
        self.the_structure.atoms()[self.acceptor_o].pos.distance_to(&n_atom.pos)
    }

    /// Calculates the planar angle between the acceptor, the donated proton and the donor atom
    #[allow(non_snake_case)]
    pub fn angle_AHD(&self) -> f64 {
        let n_atom = &self.the_structure.atoms()[self.donor_n];
        let o_atom = &self.the_structure.atoms()[self.acceptor_o];
        planar_angle3(&n_atom.pos, &self.donated_h, &o_atom.pos)
    }

    /// Calculates the planar angle between the atom preceding the acceptor, the acceptor atom and the donated proton.
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
    /// E = 0.084 { 1/r_{ON} + 1/r_{CH} - 1/r_{OH} - 1/r_{CN} } * 332 kcal/mol
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

    /// The number of residues in the structure indexed by this map.
    ///
    /// Note, that not all the residues create hydrogen bonds.
    pub fn n_residues(&self) -> usize { self.indexer.n_residues() }

    /// Iterates over all the residues in the structure indexed by this [`BackboneHBondMap`].
    pub fn residue_ids(&self) -> impl Iterator<Item=&ResidueId> { self.indexer.residue_ids() }

    pub fn h_bonds(&self)-> impl Iterator<Item = (&(usize, usize), &BackboneHBond<'a>)> {
        self.h_bonds.iter()
    }

    /// Checks if the specified residue accepts a proton in an N-turn hydrogen bond.
    ///
    /// In such a case, the hydrogen bond is formed with the residue at `which_residue + n`.
    /// Proper values of n are:
    /// - n = 3: part of a $3_{10}$ helix
    /// - n = 4: part of an $\alpha$-helix
    /// - n = 5: part of a $\pi$-helix
    ///
    /// # Arguments
    /// * `which_residue` - The index of the residue of interest.
    /// * `n` - The turn type; should be 3, 4, or 5.
    ///
    /// # Returns
    /// `true` if the specified residue forms an N-turn hydrogen bond, `false` otherwise.
    pub fn accepts_n_turn(&self, acceptor_residue: &ResidueId, n: usize) -> bool {
        let a = self.indexer.index(acceptor_residue).unwrap();
        if a+n >= self.n_residues() { return false; }
        self.h_bond_for_indexes(a + n, a).is_some()
    }

    /// Checks if the specified residue accepts a proton in an N-turn hydrogen bond.
    ///
    /// In such a case, the hydrogen bond is formed with the residue at `which_residue - n`.
    /// Proper values of n are:
    /// - n = 3: part of a $3_{10}$ helix
    /// - n = 4: part of an $\alpha$-helix
    /// - n = 5: part of a $\pi$-helix
    ///
    /// # Arguments
    /// * `which_residue` - The index of the residue of interest.
    /// * `n` - The turn type; should be 3, 4, or 5.
    ///
    /// # Returns
    /// `true` if the specified residue forms an N-turn hydrogen bond, `false` otherwise.
    pub fn donates_n_turn(&self, donor_residue: &ResidueId, n: usize) -> bool {
        let d = self.indexer.index(donor_residue).unwrap();
        if d < n { return false; }
        self.h_bond_for_indexes(d, d - n).is_some()
    }

    /// Checks if two residues form a parallel beta bridge.
    ///
    /// The order of residues is unimportant in this case.
    ///
    /// # Arguments
    /// * `i_residue` - the first residue of interest.
    /// * `j_residue` - the second residue of interest.
    ///
    /// # Returns
    /// `true` if `i_residue` and `j_residue` form a parallel beta bridge, `false` otherwise.
    pub fn is_parallel_bridge(&self, i_residue: &ResidueId, j_residue: &ResidueId) -> bool {
        let i_idx = self.indexer.index(i_residue).unwrap();
        let j_idx = self.indexer.index(j_residue).unwrap();

        if i_idx > 0 && i_idx+1 < self.n_residues() {
            if self.is_hb(j_idx,i_idx-1) && self.is_hb(i_idx+1,j_idx) { return true; }
        }

        if j_idx > 0 && j_idx + 1 < self.n_residues() {
            if self.is_hb(i_idx,j_idx-1) && self.is_hb(j_idx+1,i_idx) { return true; }
        }
        return false;
    }

    /// Checks if two residues form an antiparallel beta bridge.
    ///
    /// The order of residues is unimportant in this case.
    ///
    /// # Arguments
    /// * `i_residue` - the first residue of interest.
    /// * `j_residue` - the second residue of interest.
    ///
    /// # Returns
    /// `true` if `i_residue` and `j_residue` form an antiparallel beta bridge, `false` otherwise.
    pub fn is_antiparallel_bridge(&self, i_residue: &ResidueId, j_residue: &ResidueId) -> bool {
        let i_idx = self.indexer.index(i_residue).unwrap();
        let j_idx = self.indexer.index(j_residue).unwrap();

        if self.is_hb(j_idx, i_idx) && self.is_hb(i_idx, j_idx) { return true; }
        if i_idx > 0 && j_idx > 0 && i_idx < self.n_residues() - 1 && j_idx < self.n_residues() - 1 {
            if self.is_hb(j_idx + 1, i_idx - 1) && self.is_hb(i_idx + 1, j_idx - 1) { return true; }
        }

        return false;
    }

    /// Provides an H-bond between two residues.
    ///
    /// If the residues are not H-bonded, returns `None`
    pub fn h_bond(&self, donor_residue: &ResidueId, acceptor_residue: &ResidueId) -> Option<&BackboneHBond<'a>> {

        let d = self.indexer.index(donor_residue);
        let a = self.indexer.index(acceptor_residue);
        if let (Some(d), Some(a)) = (d, a) {
            return self.h_bonds.get(&(d, a));
        } else { return None }
    }

    /// Returns ``true`` if the specified residue ``a`` accepts a hydrogen bond from a donor ``d``
    fn is_hb(&self, d: usize, a: usize) -> bool { self.h_bond_for_indexes(d, a).is_some() }

    fn find_donors_acceptors(&mut self) {

        // ---------- Extract backbone atoms ------------
        let mut bb: Vec<BackboneResidue> = vec![BackboneResidue::new(); self.indexer.n_residues()];
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
        let mut hydrogens: Vec<Option<Vec3>> = vec![None; self.indexer.n_residues()];
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