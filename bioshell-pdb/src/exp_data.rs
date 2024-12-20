use std::fmt;
use bioshell_cif::CifData;

/// Represents experimental method types used to determine a structure of a macromolecular assembly.
///
#[derive(Debug, PartialEq, Clone)]
pub enum ExperimentalMethod {
    /// X-ray diffraction method
    XRay,
    /// Fiber diffraction method
    FiberDiffraction,
    /// Neutron diffraction method
    NeutronDiffraction,
    /// Electron crystallography
    ElectronCrystallography,
    /// Electron microscopy
    ElectronMicroscopy,
    /// Solution nuclear magnetic resonance (NMR)
    SolidStateNMR,
    /// Solution nuclear magnetic resonance (NMR)
    SolutionNMR,
    /// Solution scattering
    SolutionScattering,
    /// Infrared spectroscopy
    InfraredSpectroscopy,
    /// Theoretical model
    TheoreticalModel,
}

/// Implement the Display trait for ExperimentalMethod enum
impl fmt::Display for ExperimentalMethod {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// Converts the [`ExperimentalMethod`](ExperimentalMethod) enum variants to strings
impl ExperimentalMethod {
    pub fn to_string(&self) -> &'static str {
        match self {
            ExperimentalMethod::XRay => "X-RAY DIFFRACTION",
            ExperimentalMethod::SolutionNMR => "SOLUTION NMR",
            ExperimentalMethod::ElectronMicroscopy => "ELECTRON MICROSCOPY",
            ExperimentalMethod::NeutronDiffraction => "NEUTRON DIFFRACTION",
            ExperimentalMethod::TheoreticalModel => "THEORETICAL MODEL",
            ExperimentalMethod::FiberDiffraction => "FIBER DIFFRACTION",
            ExperimentalMethod::ElectronCrystallography => "ELECTRON CRYSTALLOGRAPHY",
            ExperimentalMethod::SolidStateNMR => "SOLID-STATE NMR",
            ExperimentalMethod::SolutionScattering => "SOLUTION SCATTERING",
            ExperimentalMethod::InfraredSpectroscopy => {"INFRARED SPECTROSCOPY"}
        }
    }

    /// Creates a [`ExperimentalMethod`](ExperimentalMethod) enums from a relevant PDB-formatted line.
    ///
    /// Note, that a single EXPDTA entry may provide more than just a single experimental method;
    /// therefore a vector of [`ExperimentalMethod`](ExperimentalMethod) enums is returned
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{ExperimentalMethod};
    /// let expdata_line = "EXPDTA    ELECTRON MICROSCOPY";
    /// let method = ExperimentalMethod::from_expdata_line(expdata_line);
    /// assert_eq!(method[0].to_string(), "ELECTRON MICROSCOPY");
    /// let expdata_line2 = "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION";
    /// let methods = ExperimentalMethod::from_expdata_line(expdata_line2);
    /// assert_eq!(methods.len(), 2);
    ///
    /// ```
    pub fn from_expdata_line(line: &str) -> Vec<ExperimentalMethod> {
        let mut methods = Vec::new();
        let extracted_line = if line.starts_with("EXPDTA") { &line[10..] } else { line };
        let expdata_entries: Vec<&str> = extracted_line.trim().split(';').map(|s| s.trim()).collect();
        for entry in expdata_entries {
            let entry_only = entry.replace("\'", "");
            let tokens: Vec<&str> = entry_only.split_whitespace().collect();
            if let Some(method) = match tokens.as_slice() {
                ["X-RAY", "DIFFRACTION"] => Some(ExperimentalMethod::XRay),
                ["FIBER", "DIFFRACTION"] => Some(ExperimentalMethod::FiberDiffraction),
                ["NEUTRON", "DIFFRACTION"] => Some(ExperimentalMethod::NeutronDiffraction),
                ["ELECTRON", "CRYSTALLOGRAPHY"] => Some(ExperimentalMethod::ElectronCrystallography),
                ["ELECTRON", "MICROSCOPY"] => Some(ExperimentalMethod::ElectronMicroscopy),
                ["SOLUTION", "NMR"] => Some(ExperimentalMethod::SolutionNMR),
                ["SOLID-STATE", "NMR"] => Some(ExperimentalMethod::SolidStateNMR),
                ["SOLUTION", "SCATTERING"] => Some(ExperimentalMethod::SolutionScattering),
                ["THEORETICAL", "MODEL"] => Some(ExperimentalMethod::TheoreticalModel),
                ["INFRARED", "SPECTROSCOPY"] => Some(ExperimentalMethod::InfraredSpectroscopy),
                _ => None,
            } {
                methods.push(method);
            }
        }
        return methods;
    }

    pub fn from_cif_data(cif_data_block: &CifData) -> Vec<ExperimentalMethod> {
        let mut output: Vec<ExperimentalMethod> = vec![];

        if let Some(methods) = cif_data_block.get_item::<String>("_exptl.method") {
            output.extend(ExperimentalMethod::from_expdata_line(&methods));
        }
        if let Some(exp_loop) = cif_data_block.get_loop("_exptl.method") {
            let idx = exp_loop.column_index("_exptl.method").unwrap();
            for row in exp_loop.rows() {
                output.extend(ExperimentalMethod::from_expdata_line(&row[idx]));
            }
        }

        return output;
    }
}
