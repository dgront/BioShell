#[allow(non_upper_case_globals)]
const cif_2fdo:  &str = include_str!("./input_files/2fdo.cif");
#[allow(non_upper_case_globals)]
const cif_2gb1:  &str = include_str!("./input_files/2gb1.cif");
#[allow(non_upper_case_globals)]
const cif_1c5n:  &str = include_str!("./input_files/1c5n.cif");

#[cfg(test)]
mod test_bb_hbond_map {
    use std::io::BufReader;
    use bioshell_interactions::BackboneHBondMap;
    use bioshell_pdb::{Deposit, PDBError, ResidueId};
    use crate::cif_2gb1;

    #[test]
    fn hbonds_2gb1() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_2gb1.as_bytes());
        let strctr = Deposit::from_cif_reader(reader)?.structure()?;
        let hbonds = BackboneHBondMap::new(&strctr);

        assert!(hbonds.is_antiparallel_bridge(&ResidueId::new("A", 16, ' '),
                                      &ResidueId::new("A", 5, ' ')));
        assert!(!hbonds.is_antiparallel_bridge(&ResidueId::new("A", 16, ' '),
                                              &ResidueId::new("A", 4, ' ')));
        assert!(hbonds.is_antiparallel_bridge(&ResidueId::new("A", 14, ' '),
                                              &ResidueId::new("A", 7, ' ')));
        assert!(hbonds.is_antiparallel_bridge(&ResidueId::new("A", 15, ' '),
                                              &ResidueId::new("A", 6, ' ')));

        assert!(hbonds.donates_n_turn(&ResidueId::new("A", 36,' '), 4));

        Ok(())
    }
}

#[cfg(test)]
mod test_dssp {
    use std::io::BufReader;
    use bioshell_interactions::{BackboneHBondMap, dssp};
    use bioshell_pdb::{Deposit, PDBError};
    use super::*;

    #[test]
    fn test_dssp() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_1c5n.as_bytes());
        let strctr = Deposit::from_cif_reader(reader)?.structure()?;
        let exp_ss = strctr.secondary("H").to_string();

        let hbonds = BackboneHBondMap::new(&strctr);
        let result = dssp(&hbonds);

        assert_eq!(result, exp_ss);

        Ok(())
    }
}