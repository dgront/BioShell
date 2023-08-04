// todo: This struct must be removed from here; this functionality must be implemented in bioshell-sequences crate (currently in bioshell-core), probably as ExtendedSequence or AnnotatedSequence type
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

    pub fn get_chain_id(&self) -> &str {
        &self.chain_id
    }

    pub fn get_sequence(&self) -> &[String] {
        &self.sequence
    }
}
