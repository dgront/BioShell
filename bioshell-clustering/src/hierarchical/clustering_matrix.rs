
pub(crate) struct HierarchicalClusteringMatrix {
    pub order: usize,
    pub sizes: Vec<usize>,
    pub matrix: Vec<Vec<f32>>,
}

impl HierarchicalClusteringMatrix {

    /// Create a new distance matrix from a given distance function.
    pub fn new<F: Fn(usize, usize) -> f32>(n_data: usize, distance_function: F) -> Self {

        let mut matrix = vec![vec![0.0_f32; n_data]; n_data];
        for i in 1..n_data {
            for j in 0..i {
                let v = distance_function(i, j);
                matrix[i][j] = v;
                matrix[j][i] = v;
            }
        }
        Self { order: n_data, sizes: vec![1; n_data], matrix }
    }

    /// Find the two elements for which the matrix holds the smallest distance value.
    ///
    /// Returns a tuple (i,j) of the indices of the elements, such that i < j.
    pub fn closest_elements(&self) -> (usize, usize) {

        let mut min_distance = f32::MAX;
        let mut min_i = 0;
        let mut min_j = 0;
        for j in 1..self.order {
            for i in 0..j {
                if self.matrix[i][j] < min_distance {
                    min_distance = self.matrix[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }
        (min_i, min_j)
    }

    /// Replace the column and row at index j with the last element in the matrix.
    ///
    /// Rows are swapped, columns are copied. This means that the last row has the new values while the column is not updated.
    pub fn replace_with_last(&mut self, j: usize) {
        self.order -= 1;
        self.matrix.swap(j, self.order);
        for i in 0..self.order {
            self.matrix[i][j] = self.matrix[i][self.order];
        }
    }

    pub fn update_distances<M>(&mut self, i: usize, j: usize, merging_rule: &M, new_index: usize)
    where
        M: Fn(usize, usize, usize, f32, f32, f32) -> f32
    {
        let row_i = &self.matrix[i];
        let row_j = &self.matrix[j];
        let size_i = self.sizes[i];
        let size_j = self.sizes[j];

        let mut result = vec![0.0_f32; self.order];
        for k in 0..self.order {
            result[k] = merging_rule(size_i, size_j, self.sizes[j], self.matrix[i][j], row_i[k], row_j[k]);
        }
        for k in 0..self.order {
            self.matrix[new_index][k] = result[k];
            self.matrix[k][new_index] = result[k];
        }
        self.matrix[new_index][new_index] = 0.0;
        self.sizes[new_index] = size_j + size_i;
    }
}

