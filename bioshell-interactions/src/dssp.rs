use bioshell_pdb::ResidueId;
use crate::{BackboneHBond, BackboneHBondMap};


/// Implements the DSSP algorithm to detect the secondary structure of a protein.
///
/// # Arguments
/// * `n_residues` - Total number of residues in the protein.
///
/// # Returns
/// A String representing the secondary structure, where each character corresponds to a residue:
/// * 'H' - Alpha-helix
/// * 'G' - 3-10 helix
/// * 'I' - Pi-helix
/// * 'E' - Beta-strand
/// * 'C' - Coil (not part of a helix or a sheet)
pub fn dssp(bb_h_bonds: &BackboneHBondMap) -> String {

    // count the residues in the structure
    let n_residues = bb_h_bonds.n_residues();

    // Initialize the secondary structure for each residue as 'C' (coil) by default
    let mut secondary_structure: Vec<char> = vec!['C'; n_residues];

    let res_ids: Vec<&ResidueId> = bb_h_bonds.residue_ids().collect();

    // Identify 3-10 helices using the i, i+3 hydrogen bond pattern
    for i in 1..n_residues.saturating_sub(3) {
        if bb_h_bonds.accepts_n_turn(res_ids[i-1], 3)&& bb_h_bonds.accepts_n_turn(res_ids[i], 3) {
            secondary_structure[i] = 'G';
            secondary_structure[i + 1] = 'G';
            secondary_structure[i + 2] = 'G';
        }
    }

    // Identify pi-helices using the i, i+5 hydrogen bond pattern
    for i in 1..n_residues.saturating_sub(5) {
        if bb_h_bonds.accepts_n_turn(res_ids[i-1], 5) && bb_h_bonds.accepts_n_turn(res_ids[i], 5) {
            secondary_structure[i] = 'I';
            secondary_structure[i + 1] = 'I';
            secondary_structure[i + 2] = 'I';
            secondary_structure[i + 3] = 'I';
            secondary_structure[i + 4] = 'I';
        }
    }

    // Identify alpha-helices using the i, i+4 hydrogen bond pattern
    for i in 1..n_residues.saturating_sub(4) {
        if bb_h_bonds.accepts_n_turn(res_ids[i-1], 4) && bb_h_bonds.accepts_n_turn(res_ids[i], 4) {
            secondary_structure[i] = 'H';
            secondary_structure[i + 1] = 'H';
            secondary_structure[i + 2] = 'H';
            secondary_structure[i + 3] = 'H';
        }
    }

    // Step 4: Identify beta-sheets using hydrogen bonds between distant residues
    for i in 0..n_residues {
        if secondary_structure[i] == 'H' || secondary_structure[i] == 'G' || secondary_structure[i] == 'I' { continue; }
        let res_i = res_ids[i];
        for j in 0..n_residues {
            if secondary_structure[j] == 'H' || secondary_structure[j] == 'G' || secondary_structure[j] == 'I' { continue; }
            let res_j = res_ids[j];
            if bb_h_bonds.is_antiparallel_bridge(res_i, res_j) {
                secondary_structure[i] = 'E';
                secondary_structure[j] = 'E';
            }
            if bb_h_bonds.is_parallel_bridge(res_i, res_j) {
                secondary_structure[i] = 'E';
                secondary_structure[j] = 'E';
            }
        }
    }

    // Convert the secondary structure to a string
    secondary_structure.into_iter().collect()
}


