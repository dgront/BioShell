use crate::{CifData, CifError, CifLoop};

/// Provides access to data stored in selected columns of a [`CifLoop`](CifLoop).
///
/// The purpose of this struct is to simplify access data stored in a loop block. It allows to iterate over
/// rows of a loop block and unpack values of selected columns at the same time.
///
/// ```
/// use std::io::BufReader;
/// use bioshell_cif::{read_cif_buffer, CifLoop, CifError, CifTable};
///
/// # fn main() -> Result<(), CifError> {
/// let cif_data = "data_coords
/// loop_
/// _atom_site_label
/// _atom_site_x
/// _atom_site_y
/// _atom_site_z
/// _atom_site_charge
/// O1 4.154 5.699 3.026 0.0
/// C2 5.630 5.087 4.246 0.0";
///
/// let data_block = &read_cif_buffer(&mut BufReader::new(cif_data.as_bytes()))?[0];
///
/// let cif_table = CifTable::new(data_block, "_atom_site_", ["_atom_site_x", "_atom_site_y", "_atom_site_z"])?;
/// for [x, y, z] in cif_table.iter() {
///     println!("{} {} {}", x, y, z);
/// }
/// # Ok(())
/// # }
/// ```
/// It allows to provide column names only partially, as long as they are unique:
/// ```
/// use std::io::BufReader;
/// use bioshell_cif::{read_cif_buffer, CifLoop, CifError, CifTable};
/// # fn main() -> Result<(), CifError> {
/// let cif_data = "data_5edw
/// loop_
/// _pdbx_audit_revision_history.ordinal
/// _pdbx_audit_revision_history.data_content_type
/// _pdbx_audit_revision_history.major_revision
/// _pdbx_audit_revision_history.minor_revision
/// _pdbx_audit_revision_history.revision_date
/// 1 'Structure model' 1 0 2016-11-02
/// 2 'Structure model' 1 1 2017-11-22
/// 3 'Structure model' 1 2 2023-09-27";
/// let data_block = &read_cif_buffer(&mut BufReader::new(cif_data.as_bytes()))?[0];
/// let rev_table = CifTable::new(data_block, "_pdbx_audit_revision_history", ["revision_date"])?;
/// for [date] in rev_table.iter() {
///     println!("{}", date);
/// }
/// # assert_eq!(rev_table.iter().count(), 3);
///
/// # Ok(())
/// # }
/// ```
///
/// ```
/// use std::io::BufReader;
/// use bioshell_cif::{CifTable, read_cif_buffer, CifError};
/// # fn main() -> Result<(), CifError> {
/// let cif_data = "data_coords
/// _atom_site_x    4.154
/// _atom_site_y    5.699
/// _atom_site_z    3.026";
/// let data_block = &read_cif_buffer(&mut BufReader::new(cif_data.as_bytes()))?[0];
/// let xyz_table = CifTable::new(data_block, "_atom_", ["_atom_site_x", "_atom_site_y", "_atom_site_z"])?;
/// for [x, y, z] in xyz_table.iter() {
///     println!("{} {} {}", x, y, z);
/// }
/// # Ok(())
/// # }
/// ```
pub struct CifTable<'a, const N: usize> {
    cif_loop: Option<&'a CifLoop>,
    column_indices: [usize; N],         // Use a constant-sized array for the indices
    data_items: Option<[String; N]>,   // Data items if we don't have a loop block
}

impl<'a, const N: usize> CifTable<'a, N> {
    pub fn new(cif_data_block: &'a CifData, selected_loop: &str, selected_columns: [&str; N]) -> Result<Self, CifError> {

        if let Some(cif_loop) = cif_data_block.get_loop(selected_loop) {
            // Find the indices of the selected column names in the original loop
            let mut column_indices = [0; N];
            for (i, &col_name) in selected_columns.iter().enumerate() {
                let pos = cif_loop.column_names.iter().position(|name| name.contains(col_name));
                if let Some(index) = pos {
                    column_indices[i] = index;
                } else {
                    return Err(CifError::MissingCifLoopKey { item_key: col_name.to_string() });
                }
            }

            return Ok(CifTable { cif_loop: Some(cif_loop), column_indices, data_items: None, });
        } else {
            let mut data_items = [const {String::new()}; N];
            for (column_index, column_name) in selected_columns.iter().enumerate() {
                let item_name = if column_name.starts_with(selected_loop) {
                    column_name.to_string()
                } else {
                    format!("{}{}", selected_loop, column_name)
                };
                if let Some(item_value) = cif_data_block.get_item(&item_name) {
                    data_items[column_index] = item_value;
                } else {
                    return Err(CifError::MissingCifLoopKey { item_key: item_name });
                }
            }

            return Ok(CifTable { cif_loop: None, column_indices: [0; N], data_items: Some(data_items), });
        }
    }

    // Method to return an iterator over the filtered rows, borrowing from the table
    pub fn iter(&'a self) -> CifTableIter<'a, 'a, N> {
        CifTableIter {
            cif_table: self,
            row_index: 0,
            internal_array: [""; N],  // Internal array allocated once
            data_items_used: false,
        }
    }
}

/// Iterator for CifTable provides a reference to an array of &str containing the requested data entries
pub struct CifTableIter<'a, 'b, const N: usize> {
    cif_table: &'b CifTable<'a, N>,  // Borrowing from CifTable
    row_index: usize,
    internal_array: [&'a str; N],    // Internal array to store references with the same lifetime as the loop data
    data_items_used: bool,
}

impl<'a, 'b, const N: usize> Iterator for CifTableIter<'a, 'b, N>  where 'b: 'a, {
    type Item = [&'a str; N];  // Return a reference to an array of &str of length N

    fn next(&mut self) -> Option<Self::Item> {

        if let Some(values) = &self.cif_table.data_items {
            if self.data_items_used { return None; }
            self.data_items_used = true;
            for (i, value) in values.iter().enumerate() {
                self.internal_array[i] = value;
            }
            return Some(self.internal_array);
        }

        let a_loop = self.cif_table.cif_loop.unwrap();
        if self.row_index >= a_loop.data_rows.len() {
            return None;  // No more rows to iterate
        }

        // Get the current row
        let current_row = &a_loop.data_rows[self.row_index];
        self.row_index += 1;

        // Populate the internal array with references to the selected columns
        for (i, &col_idx) in self.cif_table.column_indices.iter().enumerate() {
            self.internal_array[i] = &current_row[col_idx];
        }

        // Return a reference to the internal array
        Some(self.internal_array)
    }
}



