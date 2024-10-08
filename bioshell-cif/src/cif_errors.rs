use thiserror::Error;

/// Errors that may result while using the bioshell-cif crate
#[derive(Debug, Error)]
pub enum CifError {

    #[error("General I/O error occurred while reading a CIF file")]
    /// I/O error occurred while reading a PDB or a CIF file
    Io(#[from] std::io::Error),

    #[error("Data key of a loop block not found in CIF input: {item_key}")]
    /// A loop block miss a mandatory column identified by the key: `item_key`
    MissingCifLoopKey {item_key: String},

    #[error("CIF input contains more than one data block")]
    /// CIF input contains more than one data block
    ExtraDataBlock,

    #[error("CIF input must contain at least one data block!")]
    /// CIF input must contain at least one data block!
    NoDataBlock,

    #[error("Can't parse the following CIF line: {line}")]
    /// Can't parse a CIF
    InvalidCifLine {line: String },

    #[error("The following mandatory data item can't be found in a CIF input: {data_name}")]
    /// Required data item is missing in a CIF input
    MissingCifDataKey {data_name: String},

    #[error("Can't parse the following string: {item} into {type_name} type:\n{details}")]
    /// The `item` string can't be parsed into a variable of `type_name` type
    ItemParsingError {item: String, type_name: String, details: String},

    #[error("Can't insert a data name {data_name} to a loop block which already contains data values")]
    /// Can't insert a data name to a loop block which already contains data values
    MisplacedDataNameInLoop {data_name: String},

    #[error("Can't process a data item name: {data_name}")]
    /// Can't process a data item name
    DanglingDataItem {data_name: String},

    #[error("Data values found outside of a loop block:\n {breaking_line}")]
    /// Data values found outside a loop block
    DataValuesOutsideLoop { breaking_line: String },

    #[error("The number of elements in a loop row {n_rows} is greater to the number of columns {n_columns}")]
    /// The number of elements in a loop row is greater to the number of columns
    TooManyDataValues { n_rows: usize, n_columns: usize },

    #[error("Multiline data value {data_value} found their no loop or data entry open to store it")]
    /// Multiline data value found their no loop or data entry open to store it
    MultilineStringOutsideDataItem { data_value: String },

}