
#[cfg(test)]
mod test_geometric_calculations {
    use bioshell_pdb::{PdbAtom, ResidueId, Structure};
    use bioshell_pdb::calc::{dihedral_angle4, distance, planar_angle2, planar_angle3, Vec3};

    #[allow(non_upper_case_globals)]
    const pdb_2gb1:  &str = include_str!("./test_files/2gb1.pdb");

    #[test]
    fn test_distances() {
        let lines: Vec<_> = pdb_2gb1.split("\n").filter(|&l| l.starts_with("ATOM")).collect();
        let atoms: Vec<PdbAtom> = lines.iter().map(|l| PdbAtom::from_atom_line(l)).collect();
        let strctr = Structure::from_iterator("1xyz", atoms.iter().cloned());
        let ai = strctr.atom(&ResidueId::new("A", 14, ' '), " CA ").unwrap();
        let aj = strctr.atom(&ResidueId::new("A", 15, ' '), " CA ").unwrap();
        assert_eq!(ai.res_name, "GLY");
        assert_eq!(aj.res_name, "GLU");
        assert!((distance(ai, aj) - 3.8).abs() < 0.1);
    }

    #[test]
    fn calculate_planar_angles() {
        let v0 = Vec3::new(1.5, 0.0, 0.0);
        let v1 = Vec3::new(0.0, 1.5, 0.0);
        assert!(f64::abs(planar_angle2(&v0, &v1) - std::f64::consts::PI / 2.0) < 1e-5);

        // triangle 60 deg.
        let d: f64 = 1.5;
        let a = Vec3::new(-d, 0.0, 0.0);
        let b = Vec3::new(0.0, d * (3.0 as f64).sqrt(), 0.0);
        let c = Vec3::new(d, 0.0, 0.0);
        assert!(f64::abs(planar_angle3(&a, &b, &c) - 1.0472) < 1e-4);
    }

    #[test]
    fn calculate_dihedral_angles() {
        // Phi angle for TRP43 of 2gb1
        let n = Vec3::new(3.501, -0.969, -8.009);
        let ca = Vec3::new(2.365, -1.045, -7.038);
        let c = Vec3::new(1.324, -2.064, -7.504);
        let o = Vec3::new(0.970, -2.104, -8.667);
        assert!(f64::abs(dihedral_angle4(&n, &ca, &c, &o) * 180.0 / std::f64::consts::PI
            + 44.01818450297304) < 1e-4);
    }
}