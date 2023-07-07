/*
pub struct PdbHeader {
    description: String,
    date: String,
    protein_name: String,
}

impl PdbHeader {
    pub fn new(header_line: &str) -> Self {
        Self {
            description: header_line[10..50].trim().to_string(),
            date: header_line[50..59].trim().to_string(),
            protein_name: header_line[59..66].trim().to_string(),
        }
    }

    pub fn to_string(&self) -> String {
        format!("{}{}{}", self.description, self.date, self.protein_name)
    }
}
 */