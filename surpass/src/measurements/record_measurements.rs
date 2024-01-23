use std::error::Error;
use std::marker::PhantomData;
use bioshell_io::out_writer;
use crate::measurements::SystemMeasurement;
use crate::SurpassAlphaSystem;

pub struct RecordMeasurements<T, M: SystemMeasurement<T>> {
    fname: String,
    measurements: Vec<M>,
    phantom: PhantomData<T>
}

impl<T: std::fmt::Display, M: SystemMeasurement<T>> RecordMeasurements<T, M> {
    pub fn new(fname: &str, measurements: Vec<M>) -> Result<RecordMeasurements<T, M>, Box<dyn Error>> {
        let mut stream = out_writer(&fname, false);
        stream.write(b"#")?;
        for o in &measurements {
            stream.write(o.header().as_bytes())?;
        }
        stream.write(b"\n");
        Ok(RecordMeasurements { fname: fname.to_string(), measurements, phantom: Default::default() })
    }

    pub fn add_measurement(&mut self, measurement: M) -> &mut Self {
        self.measurements.push(measurement);
        return &mut self;
    }

    pub fn observe(&self, system: &SurpassAlphaSystem) -> Result<(), Box<dyn Error>> {
        let mut stream = out_writer(&self.fname, true);
        for o in &self.measurements {
            let oi = o.measure(system);
            stream.write(format!("{:}", oi).as_bytes())?;
        }
        stream.write(b"\n");
        Ok(())
    }
}