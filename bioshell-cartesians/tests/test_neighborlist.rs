mod test_cartesian_system;

#[cfg(test)]
mod neighbor_list_test {
    use bioshell_cartesians::{Coordinates, NbListRules, NeighborList, PolymerRules};
   // use bioshell_numerical::Vec3;
    use bioshell_numerical::{assert_eq_float};//, assert_eq_vec3};

    #[test]
    fn neighbor_list_initialization_test() {
        const CUT_OFF: f64 = 6.0;
        const MAX_MOVE_RANGE: f64 = 1.0;
        let buffer_thickness = MAX_MOVE_RANGE * 4.0;
        let nb_rules: Box<dyn NbListRules> = Box::new(PolymerRules {});
        let nbl: NeighborList =
            NeighborList::new(CUT_OFF, buffer_thickness, nb_rules);

        let tolerance: f64 = 0.000001;
        assert_eq_float!(CUT_OFF, nbl.get_cutoff(), tolerance);
        assert_eq_float!(buffer_thickness, nbl.get_buffer_thickness(), tolerance);
        assert_eq!(0, nbl.get_neighbors().len());
        assert_eq!(0, nbl.get_recent_pos().len());
        assert_eq!(0, nbl.get_travelled().len());
    }

    #[test]
    fn neighbor_list_update_test() {
        const CUT_OFF: f64 = 6.0;
        const MAX_MOVE_RANGE: f64 = 1.0;
        let buffer_thickness = MAX_MOVE_RANGE * 4.0;
        let nb_rules: Box<dyn NbListRules> = Box::new(PolymerRules {});
        let mut nbl: NeighborList =
            NeighborList::new(CUT_OFF, buffer_thickness, nb_rules);

        let mut coords = Coordinates::new(11);//create a coordinate of length: 3.
        coords.add_xyz(0, 1.0, 1.0, 1.0);
        coords.add_xyz(1, 2.0, 2.0, 2.0);
        coords.add_xyz(2, 3.0, 3.0, 3.0);
        coords.add_xyz(3, 4.0, 4.0, 4.0);
        coords.add_xyz(4, 5.0, 5.0, 5.0);
        coords.add_xyz(5, 6.0, 6.0, 6.0);
        coords.add_xyz(6, 6.0, 6.0, 6.0);
        coords.add_xyz(7, 7.0, 7.0, 7.0);
        coords.add_xyz(8, 8.0, 8.0, 8.0);
        coords.add_xyz(9, 9.0, 9.0, 9.0);
        coords.add_xyz(10, 10.0, 10.0, 10.0);

        nbl.update_all(&mut coords);

        assert_eq!(0, nbl.get_neighbors().len());
        assert_eq!(0, nbl.get_recent_pos().len());
        assert_eq!(0, nbl.get_travelled().len());
    }

    #[test]
    fn closest_distance_test() {

    }
}
