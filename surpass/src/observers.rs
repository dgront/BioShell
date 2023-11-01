use std::error::Error;
use std::io::Write;
use bioshell_io::out_writer;
use bioshell_pdb::calc::Vec3;
use crate::{calculate_cm, SurpassAlphaSystem};

pub struct ObserveToFile<O> {
    fname: String,
    observable: O
}

impl<O: SystemObserver> ObserveToFile<O> {
    pub fn new(fname: &str, observable: O) -> Result<ObserveToFile<O>, Box<dyn Error>> {
        let mut stream = out_writer(&fname, false);
        stream.write(observable.header().as_bytes())?;
        Ok(ObserveToFile{ fname: fname.to_string(), observable })
    }

    pub fn observe(&self, system: &SurpassAlphaSystem) -> Result<(), Box<dyn Error>> {
        let mut stream = out_writer(&self.fname, true);
        stream.write(self.observable.to_string(system).as_bytes())?;
        Ok(())
    }
}

pub trait SystemObserver {
    fn to_string(&self, system: &SurpassAlphaSystem) -> String;
    fn header(&self) -> String;
}
pub struct ObserveCM {n_chains: usize}

impl ObserveCM {
    pub fn new(system: &SurpassAlphaSystem) -> ObserveCM {
        ObserveCM{n_chains: system.count_chains()}
    }
}

impl SystemObserver for ObserveCM {
    fn to_string(&self, system: &SurpassAlphaSystem) -> String {
        let mut str: Vec<String> = vec![];
        for i_chain in 0..system.count_chains() {
            let cm_vec = calculate_cm(system, i_chain);
            str.push(format!("{:8.3} {:8.3} {:8.3}",cm_vec.x, cm_vec.y, cm_vec.z));
        }
        let mut out = str.join("  ");
        out.push('\n');

        return out;
    }

    fn header(&self) -> String {
        let mut header = String::from("#   cm-X     cm-Y    cm-Z");
        for _i in 1..self.n_chains { header.push_str("   cm-X     cm-Y    cm-Z"); }
        header.push('\n');

        return header;
    }
}

pub struct REndSquared;

impl SystemObserver for REndSquared {
    fn to_string(&self, system: &SurpassAlphaSystem) -> String {
        let mut str: Vec<String> = vec![];
        for i_chain in 0..system.count_chains() {
            let chain_atoms = system.chain_atoms(i_chain);
            let r = system.distance_squared(chain_atoms.start, chain_atoms.end - 1);
            str.push(format!("{:8.3}", r));
        }
        let mut out = str.join("  ");
        out.push('\n');

        return out;
    }

    fn header(&self) -> String { String::from("#  r-end-squared\n") }
}

pub struct RgSquared;

impl SystemObserver for RgSquared {
    fn to_string(&self, system: &SurpassAlphaSystem) -> String {
        let mut str: Vec<String> = vec![];
        for i_chain in 0..system.count_chains() {
            let chain_atoms = system.chain_atoms(i_chain);
            let cm = calculate_cm(system, i_chain);
            let (mut s, mut n, mut cc) = (0.0, 0.0, 0.0);
            let mut atom = Vec3::from_float(0.0);
            for i_atom in chain_atoms.clone() {
                n += 1.0;
                system.set_ca_to_nearest_vec3(i_atom,chain_atoms.start, &mut atom);
                cc = atom.x - cm.x;
                s += cc * cc;
                cc = atom.y - cm.y;
                s += cc * cc;
                cc = atom.z - cm.z;
                s += cc * cc;
            }
            str.push(format!("{:8.3}", s/n));
        }
        let mut out = str.join("  ");
        out.push('\n');

        return out;
    }

    fn header(&self) -> String { String::from("#  gyration-radius-squared\n") }
}