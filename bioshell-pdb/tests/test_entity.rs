#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{Deposit, Entity, EntityType, PDBError};
    use bioshell_pdb::PolymerEntityType::{PolypeptideL, DNA};
    use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager};
    use bioshell_seq::chemical::StandardResidueType::{GAP, UNK};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");
    #[allow(non_upper_case_globals)]
    const cif_5edw:  &str = include_str!("./test_files/5edw.cif");
    #[allow(non_upper_case_globals)]
    const cif_1c5n:  &str = include_str!("./test_files/1c5n.cif");
    #[allow(non_upper_case_globals)]
    const cif_2fdo: &str = include_str!("./test_files/2fdo.cif");

    #[test]
    fn test_entity_parsing() {
        let test_data = [(cif_2gb1, 1), (cif_5edw, 6)];
        for (data, num_entities) in test_data.iter() {
            let reader = BufReader::new(data.as_bytes());
            let cif_data = read_cif_buffer(reader).unwrap();
            assert_eq!(cif_data.len(), 1);
            let entities = Entity::from_cif_data(&cif_data[0]).unwrap();
            assert_eq!(entities.len(), *num_entities);
        }
    }

    #[test]
    fn check_5edw_entities() {
        let reader = BufReader::new(cif_5edw.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        // --- register residue types; this is necessary since we don't parse the CIF file here and do not extract this data from it
        let ca = ResidueType::from_attrs("CA", UNK, MonomerType::NonPolymer);
        ResidueTypeManager::get().register_residue_type(ca);
        let ttp = ResidueType::from_attrs("TTP", UNK, MonomerType::NonPolymer);
        ResidueTypeManager::get().register_residue_type(ttp);
        let entities = Entity::from_cif_data(&cif_data[0]).unwrap();

        let mut keys: Vec<&str> = entities.keys().map(|k| k.as_str()).collect();
        keys.sort();
        let expected_keys = vec!["1", "2", "3", "4", "5", "6"];
        assert_eq!(keys, expected_keys);

        let mut counts: HashMap<EntityType, usize> = [
            (EntityType::Polymer(PolypeptideL), 0),
            (EntityType::Polymer(DNA), 0),
            (EntityType::NonPolymer, 0),
            (EntityType::Water, 0)].iter().cloned().collect();
        for (_,entity) in entities.iter() {
            counts.entry(entity.entity_type()).and_modify(|e| *e += 1);
        }
        assert_eq!(counts[&EntityType::Polymer(PolypeptideL)], 1);
        assert_eq!(counts[&EntityType::Polymer(DNA)], 2);
        assert_eq!(counts[&EntityType::NonPolymer], 2);
        assert_eq!(counts[&EntityType::Water], 1);
    }

    #[test]
    fn entities_5edw_structure() -> Result<(), PDBError> {
        let cif_data = include_str!("../tests/test_files/5edw.cif");
        let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;

        assert_eq!(deposit.entities().count(), 6);

        let mut protein_entities = 0;
        let mut dna_entities = 0;
        for (id, entity) in deposit.entities() {
            if entity.entity_type() == EntityType::Polymer(PolypeptideL) {
                protein_entities += 1;
                assert_eq!(entity.chain_ids(), &["A"]);
                assert_eq!(entity.entity_monomers().len(), 341);
            }
            if entity.entity_type() == EntityType::Polymer(DNA) { dna_entities += 1; }
        }
        assert_eq!(protein_entities, 1);
        assert_eq!(dna_entities, 2);
        Ok(())
    }

    #[test]
    fn entities_2fdo_structure() -> Result<(), PDBError> {

        let deposit = Deposit::from_cif_reader(cif_2fdo.as_bytes())?;

        assert_eq!(deposit.entities().count(), 2);
        let entity = deposit.entity("1");
        assert_eq!(entity.entity_type(), EntityType::Polymer(PolypeptideL));
        let strctr = deposit.structure();

        // --- amino acids as seen in an entity
        let n_all_aa_a = entity.chain_monomers("A")?.len();
        let n_aa_a = entity.chain_monomers("A")?.iter().filter(|&rt|rt.parent_type != GAP).count();
        assert_eq!(n_all_aa_a, 94);
        assert_eq!(n_aa_a, 93); // --- one is missing

        // --- amino acids as seen in a structure
        let n_aa_a_strctr = strctr.chain_residue_ids("A").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .count();
        assert_eq!(n_aa_a_strctr, 93);

        Ok(())
    }

    #[test]
    fn test_1cn5_chain_sequences() -> Result<(), PDBError> {

        let deposit = Deposit::from_cif_reader(cif_1c5n.as_bytes())?;
        let strctr = deposit.structure();
        assert_eq!(deposit.entities().count(), 7);
        let first_entity = deposit.entity("1");
        let n_aa_l = first_entity.chain_monomers("L")?.len();
        assert_eq!(n_aa_l, 36);

        let first_entity = deposit.entity("2");
        let n_aa_h = first_entity.chain_monomers("H")?.len();
        assert_eq!(n_aa_h, 259);      // 7 are missing
        let n_gap_h = first_entity.chain_monomers("H").unwrap().iter()
            .filter(|rt| rt.parent_type==GAP).count();
        assert_eq!(n_gap_h, 7);      // 7 are missing

        let n_aa_h = strctr.chain_residue_ids("H").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking()).count();
        assert_eq!(n_aa_h, 259 - 7);   // n_AA - n_GAP

        Ok(())
    }
}


