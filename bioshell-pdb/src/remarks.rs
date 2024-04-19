use std::collections::HashMap;
use std::num::ParseFloatError;

pub struct PDBRemarks {
    remarks: HashMap<String, Vec<String>>,
}

impl PDBRemarks {

    pub fn new() -> Self {
        PDBRemarks {remarks: HashMap::new(),}
    }

    /// Adds a new remark
    pub fn add_remark(&mut self, remark_line: &str) {

        if remark_line.len() < 11 { return; }
        // Extract the remark number from the line
        let remark_number = remark_line[6..10].trim();

        // Insert the remark text into the corresponding vector in the map
        self.remarks
            .entry(remark_number.to_string())
            .or_insert_with(Vec::new)
            .push(remark_line.to_string());
    }

    /// Returns true if a given remark was loaded from a PDB file.
    ///
    /// # Example
    ///
    /// ```
    /// use bioshell_pdb::PDBRemarks;
    /// let line = "REMARK 290 CRYSTALLOGRAPHIC SYMMETRY";
    /// let mut remarks = PDBRemarks::new();
    /// remarks.add_remark(&line);
    /// assert!(remarks.has_remark("290"));
    /// ```
    pub fn has_remark(&self, remark_number: &str) -> bool {
        self.remarks.contains_key(remark_number)
    }

    /// Provides remark strings for a specific remark number
    pub fn get_remark(&self, remark_number: &str) -> Option<&Vec<String>> {
        self.remarks.get(remark_number)
    }

    pub fn resolution(&self) -> Option<f64> {
        if let Some(rem2) = self.get_remark("2") {
            for ln in rem2 {
                if ln.contains("RESOLUTION") {
                    match ln[23..30].trim().parse::<f64>() {
                        Ok(res) => {return Some(res);}
                        Err(_) => {return None;}
                    }
                }
            }
        }

        return None;
    }
}