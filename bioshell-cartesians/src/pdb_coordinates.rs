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
pub fn coordinates_to_pdb(chain: &Coordinates, i_model: i16, out_fname: &str, if_append: bool) {
    let mut out_writer = out_writer(&out_fname, if_append);

    out_writer
        .write(format!("MODEL    {i_model}\n").as_bytes())
        .ok();
    let mut i_res = 0;
    let mut prev_chain: char = ' ';
    for i in 0..chain.size() {
        let chain_id = CHAINS_ORDER
            .chars()
            .nth(chain[i].chain_id as usize)
            .unwrap();
        if chain_id != prev_chain {
            prev_chain = chain_id;
            i_res = 0;
        }
        out_writer
            .write(
                format!(
                    "ATOM   {:4}{}  ALA {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C\n",
                    i,
                    " CA ",
                    chain_id,
                    i_res,
                    chain.x(i),
                    chain.y(i),
                    chain.z(i)
                )
                .as_bytes(),
            )
            .ok();
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
    let file = match File::open(input_fname) {
        Ok(file) => file,
        Err(e) => return Err(e),
    };
    let reader = BufReader::new(file);

    let mut buff: Vec<Vec3> = vec![];

    for line in reader.lines() {
        if let Some(line_ok) = line.ok() {
            if line_ok.starts_with("ATOM  ") {
                buff.push(vec_from_atom_line(&line_ok));
            }
        }
    }

    let mut coords: Coordinates = Coordinates::new(buff.len());
    coords.set_size(buff.len());
    for i in 0..buff.len() {
        coords.set(i, buff[i].x, buff[i].y, buff[i].z);
        coords[i].chain_id = buff[i].chain_id;
    }

    return Ok(coords);
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
