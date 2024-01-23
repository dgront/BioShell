use crate::SurpassAlphaSystem;

/// Measures a property of a [`SurpassAlphaSystem`](SurpassAlphaSystem)
///
/// At each call of  the [`SystemMeasurement::measure()`](SystemMeasurement::measure()) method,
/// a property of the current conformation of a [`SurpassAlphaSystem`](SurpassAlphaSystem) is
/// evaluated and returned
pub trait SystemMeasurement<T> {
    fn measure(&self, system: &SurpassAlphaSystem) -> T;
    fn header(&self) -> String;
}