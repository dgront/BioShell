
pub struct PdbCompound {
    molecule_id: i32,
    molecule_name: String,
    chain: String,
}

impl PdbCompound {
    pub fn new(molecule_id: i32, molecule_name: &str, chain: &str) -> Self {
        Self {
            molecule_id,
            molecule_name: molecule_name.to_string(),
            chain: chain.to_string(),
        }
    }

    pub fn get_molecule_id(&self) -> i32 {
        self.molecule_id
    }

    pub fn get_molecule_name(&self) -> &str {
        &self.molecule_name
    }

    pub fn get_chain_id(&self) -> &str {
        &self.chain
    }
}