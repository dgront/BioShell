use std::io::Write;
use bioshell_io::out_writer;
use crate::{cm, SurpassAlphaSystem};

pub trait Observer {
    fn observe(&self, system: &SurpassAlphaSystem);
}
pub struct ObserveCM {
    fname: String
}

impl ObserveCM {
    pub fn new(fname: &str) -> ObserveCM {
        let mut stream = out_writer(&fname, false);
        stream.write("# cm-X cm-Y cm-Z".as_bytes());
        ObserveCM{ fname: fname.to_string() }
    }
}

impl Observer for ObserveCM {
    fn observe(&self, system: &SurpassAlphaSystem) {
        let mut stream = out_writer(&self.fname, true);
        for i_chain in 0..system.count_chains() {
            let cm_vec = cm(system, i_chain);
            stream.write(format!("{} ", &cm_vec).as_bytes());
        }
        stream.write("\n".as_bytes());
    }
}