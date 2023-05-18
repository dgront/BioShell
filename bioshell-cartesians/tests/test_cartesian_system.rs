#[cfg(test)]
mod cartesian_system_test {
    use bioshell_cartesians::{CartesianSystem, Coordinates, NbList, PolymerRules};
    use bioshell_numerical::{Rototranslation, Vec3};

    fn do_modify_vector(vec: &mut Vec3)
    {
        vec.x = 0.0;
        vec.y = 0.0;
        vec.z = 0.0;
    }

    #[test]
    fn cartesian_system_test_1() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);//create a coordinate of length: 5.
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

        let mut system_coords = system.coordinates();

        assert_eq!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0));

        let mut vec_0 = system_coords[0].clone();
        let mut vec_1 = system_coords[1].clone();
        let mut vec_2 = system_coords[2].clone();

        do_modify_vector(&mut vec_0);
        do_modify_vector(&mut vec_1);
        do_modify_vector(&mut vec_2);

        system.set_vec(0, vec_0);
        system.set_vec(1, vec_1);
        system.set_vec(2, vec_2);

        assert_eq!(system.coordinates()[0], Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(system.coordinates()[1], Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(system.coordinates()[2], Vec3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn cartesian_system_test_2() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

        assert_eq!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0));

        for i in 0..3{
            let mut vec = system.coordinates()[i].clone();
            do_modify_vector(&mut vec);
            system.set_vec(i, vec);
        }

        assert_eq!(system.coordinates()[0], Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(system.coordinates()[1], Vec3::new(0.0, 0.0, 0.0));
        assert_eq!(system.coordinates()[2], Vec3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn cartesian_system_rototranslation_test() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(5);
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 1.0,2.0,3.0);
        coords.add(2, 2.0,3.0,4.0);
        coords.add(3, 3.0,4.0,5.0);
        coords.add(4, 5.0,5.0,5.0);
        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);
        assert_eq!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(system.coordinates()[1], Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(system.coordinates()[2], Vec3::new(2.0, 3.0, 4.0));
        assert_eq!(system.coordinates()[3], Vec3::new(3.0, 4.0, 5.0));
        assert_eq!(system.coordinates()[4], Vec3::new(5.0, 5.0, 5.0));
        let angle = std::f32::consts::PI;
        let roto_tran =
            Rototranslation::around_axis(&system.coordinates()[0],
                                         &system.coordinates()[4],
                                         angle.into());
        for i in 0..3{
            let mut vec = system.coordinates()[i].clone();
            roto_tran.apply_mut(&mut vec);
            system.set_vec(i, vec);
        }
        assert_eq!(system.coordinates()[1], Vec3::new(3.0, 2.0, 1.0));
        assert_eq!(system.coordinates()[2], Vec3::new(4.0, 3.0, 2.0));
        assert_eq!(system.coordinates()[3], Vec3::new(5.0, 4.0, 3.0));
    }
}