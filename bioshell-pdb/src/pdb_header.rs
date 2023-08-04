pub struct PdbHeader {
    description: String,
    date: String,
    protein_name: String,
}

impl PdbHeader {
    pub fn new(header_line: &str) -> Self {
        Self{
            description: header_line[10..50].trim().to_string(),
            date: header_line[50..59].trim().to_string(),
            protein_name: header_line[59..66].trim().to_string(),
        }
    }

    pub fn to_string(&self) -> String {
        format!(
            "{:10}{:40}{:12}{:18}",
            "HEADER",
            self.description,
            self.date,
            self.protein_name
        )
    }

    pub fn set_description(&mut self, description: &str) {
        self.description = description.to_string();
    }

    pub fn get_description(&self) -> &str {
        &self.description
    }

    pub fn set_date(&mut self, date: &str) {
        self.date = date.to_string();
    }

    pub fn get_date(&self) -> &str {
        &self.date
    }

    pub fn set_protein_name(&mut self, protein_name: &str) {
        self.protein_name = protein_name.to_string();
    }

    pub fn get_protein_name(&self) -> &str {
        &self.protein_name
    }
}