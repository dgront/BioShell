

#[cfg(test)]
mod tests {
    use itertools::assert_equal;
    use bioshell_pdb::{assert_delta, PDBRemarks};

    #[test]
    fn test_remarks_loading() {
        let line = "REMARK 290 CRYSTALLOGRAPHIC SYMMETRY";
        let mut remarks = PDBRemarks::new();
        remarks.add_remark(&line);
        assert!(remarks.has_remark("290"));
    }

    #[test]
    fn parse_resolution() {
        let line = "REMARK   2 RESOLUTION.    1.74 ANGSTROMS.";
        let mut remarks = PDBRemarks::new();
        remarks.add_remark(&line);
        if let Some(rem2) = remarks.get_remark("2") {
            for ln in rem2 {
                assert!(ln.contains("RESOLUTION"));
            }
        }
        let res_value = remarks.resolution();
        assert!(res_value.is_some());
        assert_delta!(res_value.unwrap(), 1.74, 0.0001, "incorrect resolution value");
    }
}