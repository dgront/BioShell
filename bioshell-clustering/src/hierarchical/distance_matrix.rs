use std::fmt::Display;

use data_matrix::DataMatrix;

/// Matrix of distances between elements subjected for hierarchical clustering
pub struct DataMatrixDistance {
    datamatrix: DataMatrix
}

impl DataMatrixDistance {
    pub fn from_datamatrix(datamatrix: DataMatrix) -> Self {
        Self { datamatrix }
    }

    pub fn distance(&self, i: usize, j: usize) -> f32 {
        if i == j { return 0.0; }
        self.datamatrix.get(i, j).unwrap_or(0.0) as f32
    }

    pub fn n_elements(&self) -> usize {
        self.datamatrix.nrows()
    }

    pub fn element_id(&self, i: usize) -> &str {
        &self.datamatrix.row_labels()[i]
    }
}

impl Display for DataMatrixDistance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.n_elements() {
            for j in 0..self.n_elements() {
                write!(f, "{} {} {:.2}\t", i, j, self.distance(i, j))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}