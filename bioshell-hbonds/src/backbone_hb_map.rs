use bioshell_pdb::{PdbAtom, ResidueId, Structure};
use bioshell_pdb::calc::Vec3;

pub struct BackboneHbMap<'a> {
    the_structure: Structure,
    donated_h: Vec<PdbAtom>,
    donors_n: Vec<&'a PdbAtom>,
    accpetor_o: Vec<&'a PdbAtom>,
    accpetor_c: Vec<&'a PdbAtom>,
}

impl<'a> BackboneHbMap<'a> {
    pub fn new(strctr: Structure) -> Self {
        return BackboneHbMap { the_structure: strctr,
            donated_h: vec![], donors_n: vec![],
            accpetor_o: vec![], accpetor_c: vec![]
        };
    }

    pub fn are_h_bonded(&self, res1: &ResidueId, res2: &ResidueId) -> bool {
        // let res1_backbone = self.the_structure.

        return false;
    }

    fn find_donors_acceptors(&mut self) {

        let mut prev_c : Option<Vec3> = None;
        for ri in self.the_structure.residue_ids() {
            let the_n = self.the_structure.atom(&ri, "N");
            let the_ca = self.the_structure.atom(&ri, "CA");
            let the_c = self.the_structure.atom(&ri, "C");
            let the_o = self.the_structure.atom(&ri, "O");
            if let (Some(c), Ok(n), Ok(ca)) = (prev_c, the_n, the_ca) {
                self.donated_h.push(peptide_hydrogen(&c, &n.pos, &ca.pos));
                // self.donors_n.push(n);
            }
            if let (Ok(c), Ok(o)) = (the_c, the_o) {
                prev_c = Some(c.pos);
                self.accpetor_o.push(o);
                self.accpetor_c.push(c);
            }
        }
    }
}

pub fn peptide_hydrogen(prev_c: &Vec3, the_n: &Vec3, the_ca: &Vec3) -> PdbAtom {
    return PdbAtom::new();
}