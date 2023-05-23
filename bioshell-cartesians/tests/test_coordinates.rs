mod test_cartesian_system;

#[cfg(test)]
mod coordinates_test {
    use bioshell_cartesians::Coordinates;
    use bioshell_numerical::Vec3;

    #[test]
    fn coordinates_initialization_test() {
        let mut coords = Coordinates::new(3);//create a coordinate of length: 5.
        coords.add(0, 1.0,1.0,1.0);
        coords.add(1, 2.0,2.0,2.0);
        coords.add(2, 3.0,3.0,3.0);

        assert_eq!(coords[0], Vec3::new(1.0, 1.0, 1.0));
        assert_eq!(coords[1], Vec3::new(2.0, 2.0, 2.0));
        assert_eq!(coords[2], Vec3::new(3.0, 3.0, 3.0));
    }
}
