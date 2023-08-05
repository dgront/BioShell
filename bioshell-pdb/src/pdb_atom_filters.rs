use crate::PdbAtom;

pub trait PdbAtomPredicate {
    fn check(&self, a: &PdbAtom) -> bool;
}

pub struct ChainMatches {chain_id: String}

impl PdbAtomPredicate for ChainMatches {
    fn check(&self, a: &PdbAtom) -> bool { a.chain_id == self.chain_id }
}

pub struct IsBackbone;

impl PdbAtomPredicate for IsBackbone {
    fn check(&self, a: &PdbAtom) -> bool {
        a.name == " CA " || a.name == " C  " || a.name == " N  " || a.name == " O  "
    }
}



