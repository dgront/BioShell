use crate::PdbAtom;

/// Calculate the squared distance between two atoms
pub fn distance_squared(ai: &PdbAtom, aj:&PdbAtom) -> f64 {
    let mut d = &ai.x - &aj.x;
    let mut d2 = &d * &d;
    d = &ai.y - &aj.y;
    d2 += &d * &d;
    d = &ai.z - &aj.z;
    d2 += &d * &d;
    return d2;
}

/// Calculate the distance between two atoms
pub fn distance(ai: &PdbAtom, aj:&PdbAtom) -> f64 { distance_squared(ai, aj).sqrt() }