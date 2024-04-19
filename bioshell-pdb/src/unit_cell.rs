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
}