use std::collections::HashMap;
use std::fmt::Display;
use bioshell_io::open_file;

use csv::{ReaderBuilder, StringRecord};
use crate::errors::ClusteringError;
use crate::errors::ClusteringError::InvalidDataFormat;

/// Matrix of distances between elements subjected for hierarchical clustering
pub struct DistanceMatrix {
    n_elements: usize,
    description_to_index: HashMap<String, usize>,
    similarity_matrix: Vec<Vec<f32>>,
}

fn csv_record_error(record: &StringRecord) -> ClusteringError {
    InvalidDataFormat {
        reason: "expected three entries per line: key_i, key_j, distance".to_string(),
        data: record.iter().collect::<Vec<_>>().join(","),
    }
}

impl DistanceMatrix {

    pub fn distance(&self, i: usize, j:usize) -> f32 { self.similarity_matrix[i][j] }

    pub fn n_elements(&self) -> usize { self.n_elements }

    pub fn element_id(&self, i: usize) -> &str {
        self.description_to_index.iter().find(|(_desc, index)| *index == &i).unwrap().0
    }

    /// Reads a CSV file and populates the description_to_index map and similarity_matrix.
    pub fn from_tsv(file_path: &str) -> Result<Self, ClusteringError> {
        let reader = open_file(file_path)?;

        let mut description_to_index = HashMap::new();
        let mut similarity_matrix = Vec::new();
        let mut next_index = 0;

        let mut reader = ReaderBuilder::new()
            .has_headers(false) // Specify if the CSV has headers
            .delimiter(b'\t')
            .from_reader(reader);

        // Iterate over the CSV records
        for result in reader.records() {
            let record = result?;
            let mut tokens = record.iter();

            let key1 = tokens.next().ok_or(csv_record_error(&record))?.to_string();
            let key2 = tokens.next().ok_or(csv_record_error(&record))?.to_string();
            let val = tokens.next().ok_or(csv_record_error(&record))?.to_string();

            let sim_value: f32 = match val.parse::<f32>() {
                Ok(val) => { val }
                Err(_) => { return Err(InvalidDataFormat{ reason: "can't parse to float".to_string(), data: val }) }
            };

            // Get or insert index for key1
            let index1 = *description_to_index.entry(key1.clone()).or_insert_with(|| {
                let idx = next_index;
                next_index += 1;
                similarity_matrix.push(Vec::new()); // Add a new row for this key
                for row in &mut similarity_matrix {
                    row.resize(next_index, 0.0); // Resize all rows to match the new size
                }
                idx
            });

            // Get or insert index for key2
            let index2 = *description_to_index.entry(key2.clone()).or_insert_with(|| {
                let idx = next_index;
                next_index += 1;
                similarity_matrix.push(Vec::new()); // Add a new row for this key
                for row in &mut similarity_matrix {
                    row.resize(next_index, 0.0); // Resize all rows to match the new size
                }
                idx
            });

            // Update the similarity_matrix with the sim_value
            similarity_matrix[index1][index2] = sim_value;
            similarity_matrix[index2][index1] = sim_value; // Assume symmetric similarity
        }

        Ok(Self {
            n_elements: next_index,
            description_to_index,
            similarity_matrix,
        })
    }
}

impl Display for DistanceMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.n_elements {
            for j in 0..self.n_elements {
                write!(f, "{} {} {:.2}\t", i, j, self.similarity_matrix[i][j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}