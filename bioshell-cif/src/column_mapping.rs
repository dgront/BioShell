use crate::CifError;
/// Provides easy and efficient access to data stored in [`CifLoop`](crate::CifLoop)
///
/// Column names provided by a user are converted to integer indexes accessing respective columns.
///
/// This macro fails an assertion when any of `va` coordinates differs from the corresponding
/// coordinate of `vb` vector by more than `delta`
///
/// # Example
/// ```
/// use std::io::BufReader;
/// use bioshell_cif::{cif_columns_by_name, read_cif_buffer, CifLoop, CifError};
///
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
/// let data_blocks = read_cif_buffer(&mut BufReader::new(cif_data.as_bytes())).unwrap();
/// let first_loop = data_blocks[0].loop_blocks().next().unwrap();
/// cif_columns_by_name!(Cartesians, "_atom_site_x", "_atom_site_y", "_atom_site_z",);
/// match Cartesians::new(&first_loop) {
///     Ok(get_xyz) => {
///             let mut xyz = [""; 3];
///             for row in first_loop.rows() {
///             get_xyz.data_items(row, &mut xyz);
///         }
///     }
///     Err(e) => {panic!("{}",e)}
/// }
/// ```
#[macro_export]
macro_rules! cif_columns_by_name {
    (
        $struct_name:ident,
        $(
            $value:expr,
        )*
    ) => {
        #[derive(Debug)]
        struct $struct_name {
             values: [usize; {
                let count = 0usize $(+ { let _ = stringify!($value); 1usize })*;
                count
            }],
        }

        impl $struct_name {
            pub fn new(loop_block: &CifLoop) -> Result<$struct_name, CifError> {

                let values = [
                    $(
                        match loop_block.column_index($value) {
                            Some(num) => num,
                            None => return Err(CifError::MissingCifLoopKey{item_key: $value.to_string()})
                        },
                    )*
                ];
                Ok($struct_name { values })
            }

            pub fn data_items<'a>(&self, row: &'a Vec<String>, output: &mut [&'a str])  {

                for (output_index, &selected_index) in self.values.iter().enumerate() {
                    if let Some(s) = row.get(selected_index) {
                       output[output_index] = s;
                    }
                }
            }
        }
    };
}

