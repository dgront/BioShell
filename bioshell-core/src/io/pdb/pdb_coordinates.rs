use std::fs::File;
use std::io::{self, BufReader, BufRead, Write};

use bioshell_numerical::{Vec3};
use crate::structure::Coordinates;

use crate::utils::{out_writer};


/// Writes given [`Coordinates`](Coordinates) coordinates to a PDB file.
///
/// All atoms are written as a single poly-alanine chain and named ``" CA "``
///
/// # Arguments
/// * `chain` - coordinates to be written
/// * `i_model` - model number
/// * `out_fname` - name of the output file
///
/// # Examples
/// ```rust
/// ```
pub fn coordinates_to_pdb(chain: &Coordinates, i_model: i16, out_fname: &str, if_append: bool) {
    let mut out_writer = out_writer(&out_fname, if_append);

    out_writer.write(format!("MODEL    {i_model}\n").as_bytes()).ok();
    for i in 0..chain.size() {
        out_writer.write(format!("ATOM   {:4}{}  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C\n",
                                 i+1, " CA ", i+1, chain.x(i), chain.y(i), chain.z(i)).as_bytes()).ok();
    }
    out_writer.write(b"ENDMDL\n").ok();
}

pub fn pdb_to_coordinates(input_fname: &str) -> Result<Coordinates, io::Error> {

    let file = match File::open(input_fname) {
        Ok(file) => file,
        Err(e) => return Err(e),
    };
    let reader = BufReader::new(file);

    let mut buff: Vec<Vec3> = vec!();

    for line in reader.lines() {
        if let Some(line_ok) = line.ok() {
            if line_ok.starts_with("ATOM  ") {
                buff.push(vec_from_atom_line(&line_ok));
            }
        }
    }

    let mut coords: Coordinates = Coordinates::new(buff.len());
    for i in 0..buff.len() {
        coords.set(i, buff[i].x, buff[i].y, buff[i].z);
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

    let mut v = Vec3{x: 0.0, y:0.0, z:0.0, atom_type: 0, chain_id: 0, res_type: 0};
    v.x = atom_line[30..38].trim().parse::<f64>().unwrap();
    v.y = atom_line[38..46].trim().parse::<f64>().unwrap();
    v.z = atom_line[46..54].trim().parse::<f64>().unwrap();

    return v;
}