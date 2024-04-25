use bioshell_cif::{CifData};
use bioshell_cif::CifError::MissingCifDataKey;
use crate::{PDBError, value_or_missing_key_pdb_error};
use crate::PDBError::CifParsingError;

/// A unit cell of a crystal, containing its dimensions and angles
pub struct UnitCell {
    /// a-axis dimension
    pub a: f64,
    /// b-axis dimension
    pub b: f64,
    /// c-axis dimension
    pub c: f64,
    /// alpha angle in degrees
    pub alpha: f64,
    /// beta angle in degrees
    pub beta: f64,
    /// gamma angle in degrees
    pub gamma: f64,
    /// space group symbol
    pub space_group: String,
    /// Z value
    pub z: usize,
}

impl UnitCell {
    /// Create a new `UnitCell` construct.
    /// ## Arguments
    /// * `a` - a-axis dimension
    /// * `b` - b-axis dimension
    /// * `c` - c-axis dimension
    /// * `alpha` - alpha angle in degrees
    /// * `beta` - beta angle in degrees
    /// * `gamma` - gamma angle in degrees
    /// * `space_group` - space group symbol
    /// * `z` - Z value
    pub fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64, space_group: &str, z: usize) -> Self {
        Self {
            a, b, c, alpha, beta, gamma,
            space_group: space_group.to_string(), z,
        }
    }

    /// Creates a new UnitCell struct by parsing the CRYST1 data line
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{assert_delta, UnitCell};
    /// let line1 = "CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8";
    /// let uc = UnitCell::from_cryst1_line(line1);
    /// assert_delta!(uc.a, 52.0, 0.00001, "Incorrect unit cell dimension along a axis")
    /// ```
    pub fn from_cryst1_line(line: &str) -> UnitCell {
        let a = line[6..15].trim().parse::<f64>().unwrap();
        let b = line[15..24].trim().parse::<f64>().unwrap();
        let c = line[24..33].trim().parse::<f64>().unwrap();
        let alpha = line[33..40].trim().parse::<f64>().unwrap();
        let beta = line[40..47].trim().parse::<f64>().unwrap();
        let gamma = line[47..54].trim().parse::<f64>().unwrap();
        let space_group = &line[55..66];
        let z = line[66..70].trim().parse::<usize>().unwrap();

        return UnitCell::new(a, b, c, alpha, beta, gamma, space_group, z);
    }

    /// Creates a new UnitCell struct by extracting the necessary fields from a CIF-formatted data.
    ///
    /// Returns [`MissingCifDataKey`](MissingCifDataKey) error when at least one data item key can't
    /// be found in the provided CIF block. In other words: all items shown in the example below
    /// are mandatory for this method to return ``Ok(UnitCell)``.
    ///
    /// # Examples
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::{assert_delta, PDBError, UnitCell};
    /// let cell_data = "data_cryst_cell
    ///     _cell.length_a                         58.39
    ///     _cell.length_b                         86.70
    ///     _cell.length_c                         46.27
    ///     _cell.angle_alpha                      90.00
    ///     _cell.angle_beta                       90.00
    ///     _cell.angle_gamma                      90.00
    ///     _symmetry.space_group_name_H-M         'C 1 21 1'
    ///     _cell.Z_PDB                            1
    /// ";
    /// let data_block = &read_cif_buffer(&mut BufReader::new(cell_data.as_bytes()))[0];
    /// let uc = UnitCell::from_cif_data(data_block).unwrap();
    /// assert_delta!(uc.a, 58.39, 0.00001, "Incorrect unit cell dimension along a axis");
    ///
    /// let missing_data = "data_cryst_cell
    ///     _cell.length_b                         86.70
    ///     _cell.length_c                         46.27
    ///     _cell.angle_alpha                      90.00
    ///     _cell.angle_beta                       90.00
    /// ";
    /// let data_block = &read_cif_buffer(&mut BufReader::new(missing_data.as_bytes()))[0];
    /// let uc = UnitCell::from_cif_data(data_block);
    /// assert!(uc.is_err());
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<UnitCell, PDBError> {

        let a = value_or_missing_key_pdb_error!(cif_data, "_cell.length_a", f64);
        let b = value_or_missing_key_pdb_error!(cif_data, "_cell.length_b", f64);
        let c = value_or_missing_key_pdb_error!(cif_data, "_cell.length_c", f64);
        let alpha = value_or_missing_key_pdb_error!(cif_data, "_cell.angle_alpha", f64);
        let beta = value_or_missing_key_pdb_error!(cif_data, "_cell.angle_beta", f64);
        let gamma = value_or_missing_key_pdb_error!(cif_data, "_cell.angle_gamma", f64);
        let space_group = value_or_missing_key_pdb_error!(cif_data, "_symmetry.space_group_name_H-M", String);
        let z = value_or_missing_key_pdb_error!(cif_data, "_cell.Z_PDB", usize);

        return Ok(UnitCell::new(a, b, c, alpha, beta, gamma, &space_group, z));
    }
}

