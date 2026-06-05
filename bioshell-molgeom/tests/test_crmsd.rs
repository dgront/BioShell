
#[cfg(test)]
mod test_crmsd {
    use std::io::BufReader;
    use rand::prelude::StdRng;
    use rand::{Rng, SeedableRng};
    use bioshell_core::Vec3;
    use bioshell_molgeom::align::{crmsd, crmsd_transform};
    use bioshell_pdb::{Deposit, PDBError};
    use bioshell_pdb::pdb_atom_filters::{InvertPredicate, IsHydrogen, PdbAtomPredicate};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");

    #[test]
    fn calculate_crmsd_value() -> Result<(), PDBError> {
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

        let (rms, _rt) = crmsd_transform(&pos_a, &pos_b);
        assert!(rms < 0.15, "Expected CRMSD to be less than 0.15, got {}", rms);
        assert!(rms >0.01, "Expected CRMSD too low, got {}", rms);

        let rms_only= crmsd(&pos_a, &pos_b);
        assert!((rms-rms_only).abs() <0.0001, "CRMSD from crmsd() does not match CRMSD from crmsd_transform()");

        return Ok(());
    }

    #[test]
    fn crmsd_4_points() {
        let points_a: Vec<Vec3> = [(0, 0), (1, 1), (2, 0), (3, 1)]
            .into_iter()
            .map(|(x, y)| Vec3::new(x as f64, y as f64, 0.0))
            .collect();
        let points_b: Vec<Vec3> = [(0, 0), (1, 1), (2, 0), (3, 1)]
            .into_iter()
            .map(|(x, y)| Vec3::new(x as f64, 0.0, y as f64))
            .collect();

        let (rms, rt) = crmsd_transform(&points_a, &points_b);

        let mut err = 0.0;
        for (pa, pb) in points_a.iter().zip(points_b.iter()) {
            let mut pb_copy = rt.apply_inverse(&pb);
            pb_copy -= pa;
            err += pb_copy.length_squared();
        }
        err = (err / points_b.len() as f64).sqrt();
        assert!((err-rms).abs() < 0.00001, "CRMSD computed from transformation does not match the one computed from point differences");
    }

    #[test]
    fn crmsd_random_chain() {
        let points_a: Vec<Vec3> = [(0,  0), (1,  0), (1, -1), (2, -1), (3, -1), (3,  0), (4,  0), (4,  1), (5,  1), (6,  1)]
            .into_iter()
            .map(|(x, y)| Vec3::new(x as f64, y as f64, 0.0))
            .collect();
        let mut rng = StdRng::seed_from_u64(12345);
        let mut noise = || rng.gen_range(-0.1..=0.1);

        let points_b: Vec<Vec3> = points_a.iter().map(|p| Vec3::new(p.x + noise(), p.y + noise(), noise())).collect();

        let (rms, rt) = crmsd_transform(&points_a, &points_b);
        let mut err = 0.0;
        for (pa, pb) in points_a.iter().zip(points_b.iter()) {
            let mut pa_copy = rt.apply(&pa);
            pa_copy -= pb;
            err += pa_copy.length_squared();
        }
        err = (err / points_b.len() as f64).sqrt();
        assert!((err-rms).abs() < 0.00001, "CRMSD computed from transformation does not match the one computed from point differences");
    }
}
