#[cfg(test)]
mod cartesian_system_test {
    use bioshell_cartesians::{CartesianSystem, Coordinates, NbList, PolymerRules};
    use bioshell_numerical::Vec3;

    #[test]
    fn cartesian_system_test() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

        let mut system_coords = system.coordinates();

        assert_eq!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0));

        let   vec_0 = &system.coordinates()[0];
        let   vec_1 = &system.coordinates()[1];
        let   vec_2 = &system.coordinates()[2];

        system.add(0, vec_0.x, vec_0.y, vec_0.z);

        assert_eq!(system.coordinates()[0], Vec3::new(2.0, 2.0, 2.0));
    }

    #[test]
    fn cartesian_system_test_2() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

        assert_eq!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0));

        let vec: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..3 {
            system.set(i, 3.0, 2.0, 1.0);
            assert_eq!(system.coordinates()[i], Vec3::new(3.0, 2.0, 1.0));
            system.set_vec(i, vec);
            assert_eq!(system.coordinates()[0], Vec3::new(0.0, 0.0, 0.0));
        }
    }
}
