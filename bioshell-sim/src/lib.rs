use std::ops::Range;
use std::any::Any;

/// Defines the way how a system's energy is evaluated in a BioShell's simulation.
///
pub trait Energy<S> {
    /// Evaluates the total energy of a given system, taking into account all its components
    fn energy(&self, system: &S) -> f64;

    /// Evaluates the energy of a single component (atom, spin, residue) of a given system
    fn energy_by_pos(&self, system: &S, pos: usize) -> f64;

    /// Evaluates the energy of a contiguous range of system's components.
    fn energy_by_range(&self, system: &S, range: &Range<usize>) -> f64;

    /// Returns the name of this energy function.
    /// The returned name may be used to identify this energy, e.g. to name a column in a score table
    fn name(&self) -> String;
}

/// Defines the basic properties of a simulated system.
///
/// BioShell's approach to molecular simulations assumes each modelled system comprises a number
/// of interaction centers: spins, atoms, etc.
pub trait System: Clone {
    /// Returns the current size of the modelled system
    fn size(&self) -> usize;

    /// Sets i-th component of this system by copying its DOFs from a given `rhs` object
    fn copy_from(&mut self, i:usize, rhs: &Self);
}

/// Defines a [`System`](System) that can change its number of atoms, DOFs, etc.
pub trait ResizableSystem: System {

    /// Changes the current size of the modelled system
    fn set_size(&mut self, new_size: usize);

    /// Returns the maximum size of the modelled system
    fn capacity(&self) -> usize;
}


/// Observer takes observations of a system of a generic type `O`
pub trait Observer {
    /// The type of objects being observed byt his observer
    type O;
    /// Takes observations
    fn observe(&mut self, object: &Self::O);
    fn flush(&mut self);
    fn name(&self) -> &str;
    fn as_any(&self) -> &dyn Any;
}


/// A set of observers, that observe a system of a generic type `S`
///
///
pub struct  ObserversSet<S: 'static> {
    n_called: u32,
    observers: Vec<Box<dyn Observer<O = S>>>,
    lag_times : Vec<u32>
}

impl<S> ObserversSet<S> {

    pub fn new() -> ObserversSet<S> {
        ObserversSet {n_called: 0, observers:Vec::new(), lag_times:Vec::new() }
    }

    pub fn add_observer(&mut self, o: Box<dyn Observer<O=S>>, lag_time: u32) {
        self.observers.push(o);
        self.lag_times.push(lag_time);
    }

    pub fn get_observers(&self) -> &Vec<Box<dyn Observer<O=S>>> { &self.observers }

    pub fn observe(&mut self, object: &S) {
        for i in 0..self.observers.len() {
            if self.n_called % self.lag_times[i] == 0 {
                self.observers[i].observe(object);
            }
        }
        self.n_called += 1;
    }

    /// Call `flush()` method for all observers this sampler posses
    /// This typically writes data to streams and clears buffers
    pub fn flush_observers(&mut self) { for o in self.observers.iter_mut() { o.flush(); } }

    pub fn get_observer<T: 'static>(&self, name: &str) -> Option<&T> {
        for o in self.get_observers() {
            if name == o.name() { return o.as_any().downcast_ref::<T>(); }
        }

        return None;
    }
}