use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

use crate::Coordinates;
use bioshell_numerical::Vec3;
use bioshell_sim::{ResizableSystem, System};

use bioshell_core::utils::out_writer;

const CHAINS_ORDER: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

/// Writes given [`Coordinates`](Coordinates) coordinates to a PDB file.
///
/// All atoms are written as poly-alanine chains and named ``" CA "``; chains are indexed by
/// subsequent capital letters, followed by digits if necessary. **Note**: this method can handle
/// a system comprising no more that 37 chains.
///
/// # Arguments
/// * `chain` - coordinates to be written
/// * `i_model` - model number
/// * `out_fname` - name of the output file
/// * `if_append` - if false, this function will clear the content of a file and overwrite it with a new content
///
/// # Examples
/// ```rust
/// ```
pub fn coordinates_to_pdb(coords: &Coordinates, i_model: i16, out_fname: &str, if_append: bool) {
    let mut out_writer = out_writer(&out_fname, if_append);

    out_writer.write(format!("MODEL    {i_model}\n").as_bytes()).ok();
    let mut i_res = 0;
    let mut prev_chain: char = ' ';
    for i in 0..coords.get_size() {
        let chain_id = CHAINS_ORDER
            .chars()
            .nth(coords[i].chain_id as usize)
            .unwrap();
        if chain_id != prev_chain {
            prev_chain = chain_id;
            i_res = 0;
        }
        out_writer
            .write(format!(
                "ATOM   {:4}{}  ALA {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C\n",
                i, " CA ", chain_id, i_res, coords.get_x(i), coords.get_y(i), coords.get_z(i)).as_bytes(), ).ok();
        i_res += 1
    }
    out_writer.write(b"ENDMDL\n").ok();
}

/// Reads a content of a PDB file into a  [`Coordinates`](Coordinates) object.
///
/// This function supports only single-model files, i.e. it can't be used to load a PDB trajectory
/// with multiple conformations
///
/// # Arguments
/// * `input_fname` - name of the input file
///
/// # Examples
/// ```rust
/// ```
pub fn pdb_to_coordinates(input_fname: &str) -> Result<Coordinates, io::Error> {
    let file = match File::open(input_fname) { // Open the file
        Ok(file) => file, // File opened successfully
        Err(e) => return Err(e), // Return an error if file opening fails
    };

    let reader = BufReader::new(file); // Wrap the file in a buffered reader

    let mut buff: Vec<Vec3> = vec![]; // Create an empty vector to store coordinates

    for line in reader.lines() { // Iterate over each line in the reader
        if let Some(line_ok) = line.ok() { // Check if line is valid
            if line_ok.starts_with("ATOM  ") { // Check if line starts with "ATOM  "
                buff.push(vec_from_atom_line(&line_ok)); // Parse line and add coordinates to the vector
            }
        }
    }

    let mut coords: Coordinates = Coordinates::new(buff.len()); // Create a new Coordinates instance
    coords.set_size(buff.len()); // Set the size of the Coordinates instance

    for i in 0..buff.len() { // Iterate over the coordinates vector
        coords.set_xyz(i, buff[i].x, buff[i].y, buff[i].z); // Set x, y, z coordinates for the i-th atom
        coords[i].chain_id = buff[i].chain_id; // Set the chain ID for the i-th atom
    }

    return Ok(coords); // Return the Coordinates instance
}


/// Provides X, Y, Z from an ATOM line in the pdb format
fn xzy_from_atom_line(atom_line: &String, x: &mut f64, y: &mut f64, z: &mut f64) {
    *x = atom_line[30..38].trim().parse::<f64>().unwrap();
    *y = atom_line[38..46].trim().parse::<f64>().unwrap();
    *z = atom_line[46..54].trim().parse::<f64>().unwrap();
}

/// Fills Vec3 from an ATOM line in the pdb format
fn vec_from_atom_line(atom_line: &String) -> Vec3 {
    let mut v = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
        atom_type: 0,
        chain_id: 0,
        res_type: 0,
    };
    let chars: Vec<char> = atom_line.chars().collect();
    v.chain_id = CHAINS_ORDER.find(chars[21]).unwrap() as u16;
    xzy_from_atom_line(atom_line, &mut v.x, &mut v.y, &mut v.z);
    v.x = atom_line[30..38].trim().parse::<f64>().unwrap();
    v.y = atom_line[38..46].trim().parse::<f64>().unwrap();
    v.z = atom_line[46..54].trim().parse::<f64>().unwrap();

    return v;
}
