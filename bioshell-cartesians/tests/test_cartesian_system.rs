#[cfg(test)]
mod cartesian_system_test {
    use bioshell_cartesians::{CartesianSystem, Coordinates, NbList, NbListRules, PolymerRules};
    use bioshell_numerical::{Rototranslation, Vec3, assert_eq_vec3, assert_eq_float};

    #[test]
    fn cartesian_system_test_1() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0, 1.0, 1.0);
        coords.add(1, 2.0, 2.0, 2.0);
        coords.add(2, 3.0, 3.0, 3.0);

        let system: CartesianSystem = CartesianSystem::new(coords, nbl);

        assert_eq_vec3!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0), 0.00001);
    }

    #[test]
    fn cartesian_system_test_2() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(3);   //create a coordinate of length: 3.
        coords.add(0, 1.0, 1.0, 1.0);
        coords.add(1, 2.0, 2.0, 2.0);
        coords.add(2, 3.0, 3.0, 3.0);

        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

        assert_eq_vec3!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[1], Vec3::new(2.0, 2.0, 2.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[2], Vec3::new(3.0, 3.0, 3.0), 0.00001);

        let vec: Vec3 = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..3 {
            system.set(i, 3.0, 2.0, 1.0);
            assert_eq_vec3!(system.coordinates()[i], Vec3::new(3.0, 2.0, 1.0), 0.00001);
            system.set_vec(i, vec);
            assert_eq_vec3!(system.coordinates()[0], Vec3::new(0.0, 0.0, 0.0), 0.00001);
        }
    }

    #[test]
    fn cartesian_system_rototranslation_test() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(5);
        coords.add(0, 1.0, 1.0, 1.0);
        coords.add(1, 1.0, 2.0, 3.0);
        coords.add(2, 2.0, 3.0, 4.0);
        coords.add(3, 3.0, 4.0, 5.0);
        coords.add(4, 5.0, 5.0, 5.0);
        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);
        assert_eq_vec3!(system.coordinates()[0], Vec3::new(1.0, 1.0, 1.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[1], Vec3::new(1.0, 2.0, 3.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[2], Vec3::new(2.0, 3.0, 4.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[3], Vec3::new(3.0, 4.0, 5.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[4], Vec3::new(5.0, 5.0, 5.0), 0.00001);
        let angle = std::f32::consts::PI;
        let roto_tran =
            Rototranslation::around_axis(&system.coordinates()[0],
                                         &system.coordinates()[4],
                                         angle.into());
        for i in 1..4 {
            let mut tmp = system.coordinates()[i].clone();
            roto_tran.apply_mut(&mut tmp);
            system.set_vec(i, tmp);
        }
        assert_eq_vec3!(system.coordinates()[1], Vec3::new(3.0, 2.0, 1.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[2], Vec3::new(4.0, 3.0, 2.0), 0.00001);
        assert_eq_vec3!(system.coordinates()[3], Vec3::new(5.0, 4.0, 3.0), 0.00001);
    }

    #[test]
    fn cartesian_system_periodic_image_test() {
        const E_TO: f64 = 6.0;
        let buffer_thickness = 4.0;
        let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
        let mut coords = Coordinates::new(5);
        let mut coords = Coordinates::new(5);
        coords.add(0, 1.0, 1.0, 1.0);
        coords.add(1, 1.0, 2.0, 3.0);
        coords.add(2, 2.0, 3.0, 4.0);
        coords.add(3, 3.0, 4.0, 5.0);
        coords.add(4, 5.0, 5.0, 5.0);
        let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);
        system.set_box_len(10.0);

        // Set the coordinates of a reference point and a vector
        let reference = Vec3::new(5.0, 5.0, 5.0);
        let vector = Vec3::new(11.0, -11.0, 0.0);

        // Calculate the periodic image of the vector
        let image = system.get_periodic_image(&reference, &vector);

        // Assert that the image coordinates are correctly calculated
        assert_eq!(image, Vec3::new(1.0, -1.0, 0.0));
    }
}
