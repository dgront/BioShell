use bioshell_seq::chemical::ResidueType;
use crate::{Deposit, EntityType};

const LIGAND_NOT_RELEVANT: [&str; 11] = ["HOH", "GOL", "CA", "NA", "MG", "SO4", "PO4", "MES", "CL", "K", "ZN"];

/// Lists non-trivial ligands found in a deposit
///
/// The function omits trivial molecules such as water, ions, MES and other buffer agents.
pub fn list_ligands_in_deposit(deposit: &Deposit) -> Vec<ResidueType> {
    let mut ret: Vec<ResidueType> = vec![];

    fn res_type_ok(rt: &ResidueType) -> bool {
        for rti in LIGAND_NOT_RELEVANT { if rti == rt.code3 { return false; } }
        return true;
    }

    for (_id, entity) in deposit.entities() {
        if entity.entity_type() != EntityType::NonPolymer { continue }
        for rt in entity.entity_monomers() {
            if res_type_ok(rt) {
                if ! ret.contains(rt) { ret.push(rt.clone()) };
            }
        }
    }

    return ret;
}