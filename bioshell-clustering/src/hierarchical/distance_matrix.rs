use std::fmt::Display;

use data_matrix::DataMatrix;

/// Matrix of distances between elements subjected for hierarchical clustering.
///
/// [`DataMatrixDistance`] is a wrapper around [`DataMatrix`] that
/// provides distance values loaded from a ``.tsv`` file.
///
/// # Examples
/// ```
/// use data_matrix::DataMatrixBuilder;
/// use bioshell_clustering::hierarchical::DataMatrixDistance;
/// # use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// // --- load a DataMatrix from a TSV file ---
/// let dm = DataMatrixBuilder::new()
///     .symmetric(true)
///     .default_value(0.0)
///     .from_file("tests/test_files/cities.tsv").unwrap();
///
/// // --- create a DataMatrixDistance ---
/// let dmatrix = DataMatrixDistance::from_datamatrix(dm);
/// // --- check its properties ---
/// assert_eq!(dmatrix.n_elements(), 15);
/// assert_eq!(dmatrix.distance(0, 1), 1149.36);
/// # assert_eq!(dmatrix.distance(1, 0), 1149.36);
/// # assert_eq!(dmatrix.distance(0, 0), 0.0);
/// assert_eq!(dmatrix.element_id(0), "Tokyo");
/// # assert_eq!(dmatrix.element_id(1), "Seoul");
/// # Ok(())
/// # }
/// ```
///
/// [`hierarchical_clustering()`](crate::hierarchical::hierarchical_clustering) function requires
/// a callable to compute a distance between two elements. This can be created as:
/// ```
/// # use std::error::Error;
/// # use data_matrix::DataMatrixBuilder;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// # use bioshell_clustering::hierarchical::DataMatrixDistance;
/// # let data = vec![0.0];
/// # let matrix = DataMatrixBuilder::new().from_data(&data).unwrap();
/// # let dmatrix = DataMatrixDistance::from_datamatrix(matrix);
/// let distance_fn = |i: usize, j: usize| dmatrix.distance(i, j);
/// # Ok(())
/// # }
/// ```
pub struct DataMatrixDistance {
    datamatrix: DataMatrix
}

impl DataMatrixDistance {
    /// Consumes a [`DataMatrix`] and creates a [`DataMatrixDistance`] wrapper around it.
    /// ```
    /// # use std::error::Error;
    /// # use data_matrix::DataMatrixBuilder;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// # use bioshell_clustering::hierarchical::DataMatrixDistance;
    /// # let data = vec![0.0];
    /// let matrix = DataMatrixBuilder::new().from_data(&data).unwrap();
    /// let dmatrix = DataMatrixDistance::from_datamatrix(matrix);
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_datamatrix(datamatrix: DataMatrix) -> Self {
        Self { datamatrix }
    }

    /// Returns the distance between two elements identified by their indices.
    pub fn distance(&self, i: usize, j: usize) -> f32 {
        if i == j { return 0.0; }
        self.datamatrix.get(i, j).unwrap_or(0.0) as f32
    }

    /// Returns the number of elements in the distance matrix.
    pub fn n_elements(&self) -> usize {
        self.datamatrix.nrows()
    }

    /// Returns the label of the element identified by its index.
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