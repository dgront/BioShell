macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}

const ile_pdb: [&str;8] = ["ATOM      1  N   ILE A   1       0.000   0.000   0.000  1.00  0.00           N",
    "ATOM      2  CA  ILE A   1       1.458   0.000   0.000  1.00  0.00           C",
    "ATOM      3  C   ILE A   1       2.009   1.420   0.000  1.00  0.00           C",
    "ATOM      4  O   ILE A   1       1.383   2.339  -0.529  1.00  0.00           O",
    "ATOM      5  CB  ILE A   1       2.007  -0.764  -1.218  1.00  0.00           C",
    "ATOM      6  CG1 ILE A   1       0.857  -1.301  -2.075  1.00  0.00           C",
    "ATOM      7  CG2 ILE A   1       2.916  -1.899  -0.770  1.00  0.00           C",
    "ATOM      8  CD1 ILE A   1      -0.514  -0.964  -1.536  1.00  0.00           C", ];

#[cfg(test)]
mod nerf_test {
    use bioshell_pdb::calc::{Vec3, place_atom, place_chain, dihedral_angle4, planar_angle3};
    use bioshell_pdb::PdbAtom;
    use crate::ile_pdb;

    #[test]
    fn build_stub() {
        let mut chain = vec![Vec3::default(); 4];
        let epsilon= 0.00000001;
        let right_angle = 90.0_f64.to_radians();
        let r = [0.0, 2.0, 3.0, 1.0];
        let a = [0.0, 0.0, right_angle, right_angle];
        let t = [0.0, 0.0, 0.0, right_angle];
        place_chain(&r, &a, &t, &mut chain);
        assert_delta!(chain[1].x, 2.0, epsilon);
        assert_delta!(chain[1].y, 0.0, epsilon);
        assert_delta!(chain[2].x, 2.0, epsilon);
        assert_delta!(chain[2].y, 3.0, epsilon);
        assert_delta!(chain[3].x, 2.0, epsilon);
        assert_delta!(chain[3].y, 3.0, epsilon);
        assert_delta!(chain[3].z, 1.0, epsilon);

        let t3 = dihedral_angle4(&chain[0], &chain[1], &chain[2], &chain[3]);
        assert_delta!(t[3], t3, epsilon);

        let [b, c, d] = &mut chain[1..4] else { unreachable!(); };
        let mut e = Vec3::default();
        place_atom(b, c, d, 1.0, right_angle, right_angle, &mut e);
        assert_delta!(e.x, 3.0, epsilon);
        assert_delta!(e.y, 3.0, epsilon);
        assert_delta!(e.z, 1.0, epsilon);
    }

    #[test]
    fn build_aminoacid() {

        let pdb_atoms: Vec<PdbAtom> = ile_pdb.iter().map(|&s| PdbAtom::from_atom_line(s)).collect();
        let atoms: Vec<Vec3> = pdb_atoms.iter().map(|a|a.pos.clone()).collect();

        let atom_names = [" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2"];
        let ile_topo = [[0, 0, 0, 0], [0, 1, 0, 0], [0, 1, 2, 0],
            [0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 4, 5], [1, 4, 5, 6], [1, 4, 5, 7]];

        let mut r = vec![0.0, atoms[0].distance_to(&atoms[1]), atoms[1].distance_to(&atoms[2])];
        let mut a = vec![0.0, 0.0, planar_angle3(&atoms[0],&atoms[1], &atoms[2])];
        let mut t = vec![0.0, 0.0, 0.0];
        for iatom in 3..8 {
            let [i, j, k, l] = ile_topo[iatom];
            r.push(atoms[l].distance_to(&atoms[k]));
            a.push(planar_angle3(&atoms[j],&atoms[k], &atoms[l]));
            t.push(dihedral_angle4(&atoms[i],&atoms[j], &atoms[k], &atoms[l]));
        }

        println!("{:?}", r);
        println!("{:?}", a);
        println!("{:?}", t);

        // allocate the full chain and place the first atom in its correct position
        let mut chain = vec![Vec3::default(); 8];
        chain[0].set(&atoms[0]);
        // rebuild the stub
        place_chain(&r[0..4],&a[0..4],&t[0..4],&mut chain[0..4]);
        let mut v = Vec3::default();
        for iatom in 4..8 {
            let [i, j, k, l] = ile_topo[iatom];
            place_atom(&chain[i], &chain[j], &chain[k], r[iatom], a[iatom], t[iatom], &mut v);
            chain[l].set(&v);
            println!("{}", v);
        }
            // let epsilon= 0.00000001;
        for iatom in 0..8 {
            println!("ATOM   {:4} {} ILE {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C",
                    iatom, &atom_names[iatom], "A", 1, &chain[iatom].x, &chain[iatom].y, &chain[iatom].z);
        }

    }
}