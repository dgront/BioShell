pub struct PdbSequenceOfResidue {
    chain_id: String,
    sequence: Vec<String>,
}

impl PdbSequenceOfResidue {
    pub fn new(chain_id: &str, sequence: Vec<&str>) -> Self {
        Self {
            chain_id: chain_id.to_string(),
            sequence: sequence.iter().map(|s| s.to_string()).collect(),
        }
    }

    pub fn chain_id(&self) -> &str {
        &self.chain_id
    }

    pub fn sequence(&self) -> &[String] {
        &self.sequence
    }
}