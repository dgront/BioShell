
#[cfg(test)]
mod test_crmsd {
    use std::io::BufReader;
    use rand::prelude::StdRng;
    use rand::{Rng, SeedableRng};
    use bioshell_molgeom::align::crmsd;
    use bioshell_pdb::{Deposit, PDBError};
    use bioshell_pdb::pdb_atom_filters::{InvertPredicate, IsHydrogen, PdbAtomPredicate};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");

    #[test]
    fn test_distances() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_2gb1.as_bytes());
        let deposit = Deposit::from_cif_reader(reader)?;
        let strctr = deposit.structure()?;
        let skip_h = InvertPredicate::new(IsHydrogen);
        let pos_a = strctr.atoms().iter().filter(|a| skip_h.check(a)).map(|a| a.pos.clone()).collect::<Vec<_>>();
        let mut pos_b = pos_a.iter().map(|a| a.clone()).collect::<Vec<_>>();

        let mut rng = StdRng::seed_from_u64(12345);
        let mut noise = || rng.gen_range(-0.1..=0.1);
        for p in &mut pos_b {
            p.x += noise();
            p.y += noise();
            p.z += noise();
        }

        let rms = crmsd(&pos_a, &pos_b);
        assert!(rms < 0.15, "Expected CRMSD to be less than 0.15, got {}", rms);
        assert!(rms >0.01, "Expected CRMSD too low, got {}", rms);

        return Ok(());
    }
}
