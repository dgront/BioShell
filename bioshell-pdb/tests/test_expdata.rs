

#[cfg(test)]
mod tests {
    use bioshell_pdb::{ExperimentalMethod};

    #[test]
    fn test_expdata_parsing() {
        let expdata_line = "EXPDTA    ELECTRON MICROSCOPY";
        let methods = ExperimentalMethod::from_expdata_line(expdata_line);
        assert_eq!(methods[0].to_string(), "ELECTRON MICROSCOPY");
        let expdata_line2 = "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION";
        let methods = ExperimentalMethod::from_expdata_line(expdata_line2);
        assert_eq!(methods.len(), 2);
    }
}