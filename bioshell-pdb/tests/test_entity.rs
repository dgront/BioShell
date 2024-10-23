#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{Entity, EntityType, load_cif_reader, PDBError, PolymerEntity};
    use bioshell_pdb::PolymerEntityType::{PolypeptideL, DNA};
    use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager};
    use bioshell_seq::chemical::StandardResidueType::UNK;

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
    fn load_polymer_entities() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_5edw.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let entities = PolymerEntity::from_cif_data(&cif_data[0])?;
        assert_eq!(entities.len(), 3);
        let result = vec![(341, PolypeptideL), (19, DNA), (13, DNA)];
        for entity in entities {
            let id = entity.entity_id().parse::<usize>().unwrap();
            assert_eq!(result[id-1].0, entity.monomer_sequence().len());
            assert_eq!(&result[id-1].1, entity.entity_type());
        }
        Ok(())
    }

    #[test]
    fn entities_5edw_structure() -> Result<(), PDBError> {
        let cif_data = include_str!("../tests/test_files/5edw.cif");
        let strctr = load_cif_reader(cif_data.as_bytes())?;

        assert_eq!(strctr.entities().count(), 6);

        let mut protein_entities = 0;
        let mut dna_entities = 0;
        for (id, entity) in strctr.entities() {
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

        let strctr = load_cif_reader(cif_2fdo.as_bytes())?;

        assert_eq!(strctr.entities().count(), 2);
        let entity = strctr.entity("1");
        assert_eq!(entity.entity_type(), EntityType::Polymer(PolypeptideL));

        // --- full chains, including water molecules
        assert_eq!(strctr.chain_residue_ids("B").len(), 110);
        assert_eq!(strctr.chain_residue_ids("A").len(), 104);

        // --- amino acids only
        let n_aa = strctr.chain_residue_ids("B").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .count();
        assert_eq!(n_aa, 94);
        let n_aa = strctr.chain_residue_ids("A").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .count();
        assert_eq!(n_aa, 93);
        Ok(())
    }

    #[test]
    fn test_1cn5_chain_sequences() -> Result<(), PDBError> {

        let strctr = load_cif_reader(cif_1c5n.as_bytes())?;
        assert_eq!(strctr.entities().count(), 7);
        let atom_seq: Vec<_> = strctr.chain_residue_ids("L").iter()
            .filter(|ri|strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .map(|rt| rt.clone())
            .collect();
        let entity_seq = strctr.entity("1").entity_monomers();
        let entity_chain_seq = strctr.entity("1").chain_monomers("L").unwrap();

        assert_eq!(entity_chain_seq.len(), entity_chain_seq.len());
        let entity_aa_seq: Vec<&ResidueType> = entity_chain_seq.iter().filter(|rt| rt.chem_compound_type.is_peptide_linking()).collect();
        assert_eq!(entity_aa_seq.len(), atom_seq.len());
        Ok(())
    }
}


