use std::collections::{BTreeSet, HashMap};
use crate::{PdbAtom, PDBError, Structure};
use crate::pdb_parsing_error::PdbConversionImpossibleReason;

/// Allowed PDB chain identifiers (62 total): A–Z, a–z, 0–9
static PDB_CHAIN_ID_ALPHABET: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

fn is_pdb_chain_id(s: &str) -> bool {
    s.len() == 1 && PDB_CHAIN_ID_ALPHABET.contains(s)
}

fn is_pdb_atom_serial_ok(serial: i32) -> bool {
    (1..=99_999).contains(&serial)
}

// PDB atom name field is 4 columns. We do not attempt truncation.
fn is_pdb_atom_name_ok(name: &str) -> bool {
    name.chars().count() <= 4 && !name.is_empty()
}

// Residue name is typically 3 columns in PDB.
fn is_pdb_res_name_ok(name: &str) -> bool {
    name.chars().count() <= 3 && !name.is_empty()
}

// Element symbol in PDB is 1–2 columns.
fn is_pdb_element_ok(el: &str) -> bool {
    let n = el.chars().count();
    (1..=2).contains(&n)
}

// PDB resSeq is 4 columns. While spec allows negative in theory depending on parser,
// many tools accept -999..9999. Adjust if you want stricter (e.g., 1..=9999).
fn is_pdb_res_seq_ok(res_seq: i32) -> bool {
    (-999..=9999).contains(&res_seq)
}

/// Creates a new [`Structure`] compatible with classic PDB file format.
///
/// The function clones all the atoms of the argument and odifies atom numbers, residue numbers and chain IDs where necessary.
/// The returned Structure can be safely saved in the PDB format.
pub fn make_pdb_compatible(strctr: &Structure) -> Result<Structure, PDBError> {
    let atoms = strctr.atoms();

    // ---------- Early impossibility checks ----------
    if atoms.len() > 99_999 {
        return Err(PDBError::PdbConversionNotPossible {
            reason: PdbConversionImpossibleReason::TooManyAtoms { atoms: atoms.len(), max: 62 },
        });
    }

    // Collect chain IDs as they appear (treat empty as a distinct chain label for mapping purposes).
    let mut chain_set: BTreeSet<String> = BTreeSet::new();
    for a in atoms.iter() {
        chain_set.insert(a.chain_id.clone());
    }
    if chain_set.len() > 62 {
        return Err(PDBError::PdbConversionNotPossible {
            reason: PdbConversionImpossibleReason::TooManyChains { chains: chain_set.len(), max: 62 },
        });
    }

    // ---------- Determine whether any modifications are needed ----------
    let mut needs_chain_remap = false;
    let mut needs_serial_renumber = false;
    let mut needs_res_seq_renumber = false;
    let mut needs_null_char_cleanup = false;

    // serial uniqueness detection
    let mut seen_serials: BTreeSet<i32> = BTreeSet::new();
    for a in atoms.iter() {
        // Hard validations we choose not to “fix” by truncation:
        if !is_pdb_atom_name_ok(&a.name) {
            return Err(PDBError::PdbConversionNotPossible {
                reason: PdbConversionImpossibleReason::AtomNameTooLong { atom_name: a.name.clone(), max: 4},
            });
        }
        if !is_pdb_res_name_ok(&a.res_name) {
            return Err(PDBError::PdbConversionNotPossible {
                reason: PdbConversionImpossibleReason::ResidueNameTooLong { res_name: a.res_name.clone(), max: 3}
            });
        }

        if let Some(el) = a.element.as_ref() {
            if !is_pdb_element_ok(el) {
                return Err(PDBError::PdbConversionNotPossible {
                    reason: PdbConversionImpossibleReason::ElementTooLong { element: el.clone(), max: 2},
                });
            }
        }

        // Chain IDs
        if !is_pdb_chain_id(&a.chain_id) { needs_chain_remap = true; }

        // Atom serials
        if !is_pdb_atom_serial_ok(a.serial) { needs_serial_renumber = true; }
        if !seen_serials.insert(a.serial) { needs_serial_renumber = true; }

        // Residue sequence numbers
        if !is_pdb_res_seq_ok(a.res_seq) { needs_res_seq_renumber = true; }

        // '\0' cleanups (common when coming from C-ish sources)
        if a.alt_loc == '\0' || a.i_code == '\0' { needs_null_char_cleanup = true; }
    }

    // If nothing needs changing, return an identical structure (new allocation).
    if !(needs_chain_remap || needs_serial_renumber || needs_res_seq_renumber || needs_null_char_cleanup) {
        // This assumes `PdbAtom: Clone`. If not, replace with a manual clone of fields.
        let cloned_atoms = atoms.clone();
        // Assumes you can get an id_code from Structure; adapt to your API.
        // For example: strctr.id_code()
        return Ok(Structure::from_atoms(&strctr.id_code, cloned_atoms));
    }

    // ---------- Build chain ID remapping if needed ----------
    if needs_chain_remap && chain_set.len() > 62 {
        return Err(PDBError::PdbConversionNotPossible {
            reason: PdbConversionImpossibleReason::TooManyChains { chains: chain_set.len(), max: 62}
        });
    }
    let chain_map: HashMap<String, String> = if needs_chain_remap {
        chain_set.iter()
            .zip(PDB_CHAIN_ID_ALPHABET.chars())
            .map(|(old, c)| (old.clone(), c.to_string()))
            .collect()
    } else {
        HashMap::new()
    };

    // ---------- Residue renumbering (only if needed) ----------
    //
    // We renumber per chain, preserving the order of residues as they appear in `atoms`.
    // Residue identity is tracked by (chain_id_after_mapping, entity_id, old_res_seq, i_code).
    let mut res_seq_map: HashMap<(String, String, i32, char), i32> = HashMap::new();
    let mut next_res_seq_by_chain: HashMap<String, i32> = HashMap::new();

    // ---------- Atom serial renumbering (only if needed) ----------
    let mut next_serial: i32 = 1;

    // ---------- Produce new atoms ----------
    let mut new_atoms: Vec<PdbAtom> = Vec::with_capacity(atoms.len());

    for a in atoms.iter() {
        // clone then modify as needed
        let mut na = a.clone();

        // chain_id remap
        if needs_chain_remap {
            let mapped = chain_map.get(&a.chain_id).expect("can't map chain_id even though this should work");
            na.chain_id = mapped.clone();
        }

        // clean '\0' chars
        if needs_null_char_cleanup {
            if na.alt_loc == '\0' {
                na.alt_loc = ' ';
            }
            if na.i_code == '\0' {
                na.i_code = ' ';
            }
        }

        // residue renumber
        if needs_res_seq_renumber {
            let chain_key = na.chain_id.clone();
            let key = (na.chain_id.clone(), na.entity_id.clone(), a.res_seq, a.i_code);

            let new_res_seq = if let Some(v) = res_seq_map.get(&key) {
                *v
            } else {
                let nxt = next_res_seq_by_chain.entry(chain_key).or_insert(1);
                let assigned = *nxt;
                *nxt += 1;
                res_seq_map.insert(key, assigned);
                assigned
            };
            na.res_seq = new_res_seq;
        }

        // serial renumber
        if needs_serial_renumber {
            na.serial = next_serial;
            next_serial += 1;
        }

        new_atoms.push(na);
    }

    Ok(Structure::from_atoms(&strctr.id_code, new_atoms))
}
