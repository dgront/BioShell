use std::fs;
use std::path::{MAIN_SEPARATOR};

fn file_exists(path: &str) -> bool {
    fs::metadata(path).is_ok()
}

pub(crate) fn find_deposit_file_name(pdb_code: &str,
        pdb_path: &str, prefixes: &[&str], suffixes: &[&str]) -> Result<String, std::io::Error> {

    // Check if the PDB code is actually the file we are looking for
    if file_exists(pdb_code) {
        return Ok(pdb_code.to_string());
    }

    let mut path = pdb_path.to_string();
    if !path.ends_with(MAIN_SEPARATOR) {
        path.push(MAIN_SEPARATOR);
    }

    // Check if the PDB code is the file located at the given path
    let full_path = format!("{}{}", path, pdb_code);
    if file_exists(&full_path) {
        return Ok(full_path);
    }

    let code_lo = pdb_code.to_lowercase();
    let code_up = pdb_code.to_uppercase();
    let subdir = format!("{}{}", &code_lo[1..3], MAIN_SEPARATOR);

    let mut tested = Vec::new();
    for p in prefixes {
        for s in suffixes {
            let name_lo = format!("{}{}{}{}", path, p, code_lo, s);
            if file_exists(&name_lo) {
                return Ok(name_lo);
            } else {
                tested.push(name_lo);
            }
            let name_up = format!("{}{}{}{}", path, p, code_up, s);
            if file_exists(&name_up) {
                return Ok(name_up);
            } else {
                tested.push(name_up);
            }
            let name_subdir_lo = format!("{}{}{}{}{}", path, subdir, p, code_lo, s);
            if file_exists(&name_subdir_lo) {
                return Ok(name_subdir_lo);
            } else {
                tested.push(name_subdir_lo);
            }
            let name_subdir_up = format!("{}{}{}{}{}", path, subdir, p, code_up, s);
            if file_exists(&name_subdir_up) {
                return Ok(name_subdir_up);
            } else {
                tested.push(name_subdir_up);
            }
        }
    }
    return Err(std::io::Error::new(std::io::ErrorKind::NotFound,
        format!("Could not find a  file for code '{}'. Tried: {}",
                pdb_code, tested.join(", "))));
}

mod tests_utilities {
    #[allow(unused_imports)]
    use super::*;
    #[test]
    fn test_find_file_name() {
        let pdb_code = "2gb1";
        let pdb_path = "./tests/test_files/";
        let prefixes: [&str; 4] = ["pdb", "PDB", "pdb", ""];
        let suffixes: [&str; 7] = [".ent", ".ent.gz", ".gz", ".pdb", ".PDB", ".pdb.gz", ""];

        let file_name = find_deposit_file_name(pdb_code, pdb_path, &prefixes, &suffixes);
        assert!(file_name.is_ok());
        assert_eq!(file_name.unwrap(), "./tests/test_files/2gb1.pdb");
    }
}