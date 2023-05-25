mod test_cartesian_system;

#[cfg(test)]
mod coordinates_test {
    use bioshell_cartesians::Coordinates;
    use bioshell_numerical::Vec3;
    use bioshell_numerical::{assert_eq_float, assert_eq_vec3};

    #[test]
    fn coordinates_initialization_test() {
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        assert_eq_vec3!(coords[0], Vec3::new(1.0, 1.0, 1.0), 1e-5);
        assert_eq_vec3!(coords[1], Vec3::new(2.0, 2.0, 2.0), 1e-5);
        assert_eq_vec3!(coords[2], Vec3::new(3.0, 3.0, 3.0), 1e-5);
    }

    #[test]
    fn clone_closest_image() {
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0,1.0,0.0);
        coords.add(1, 3.0,2.0,1.0);
        coords.add(2, 5.0,2.0,0.0);
        coords.set_box_len(5.0);

        let img = coords.clone_closest_image(0, 2);
        assert_eq_vec3!(img, Vec3::new(0.0, 2.0, 0.0), 1e-5);

        coords.set(0, 0.0,1.0,1.0);
        coords.set(1, 6.0,2.0,3.0);
        coords.set(2, 0.0,2.0,5.0);
        let img = coords.clone_closest_image(0, 2);
        assert_eq_vec3!(img, Vec3::new(0.0, 2.0, 0.0), 1e-5);
    }

    #[test]
    fn closest_distance() {
        let mut coords = Coordinates::new(3);//create a coordinate of length: 3.
        coords.add(0, 1.0,1.0,0.0);
        coords.add(1, 3.0,2.0,1.0);
        coords.add(2, 5.0,2.0,0.0);
        coords.set_box_len(5.0);

        let d2 = coords.closest_distance_square(0, 2);
        assert_eq_float!(d2, 2.0, 1e-5);

        let d2 = coords.closest_distance_square_to_vec(0, &Vec3::new(5.0, 2.0, 0.0));
        assert_eq_float!(d2, 2.0, 1e-5);
    }
}
