pub struct PdbTitle {
    text: String,
}

impl PdbTitle {
    pub fn new(line: &str) -> Self {
        Self {
            text: line[10..73].trim().to_string(),
        }
    }

    pub fn to_string(&self) -> String {
        self.text.clone()
    }
}