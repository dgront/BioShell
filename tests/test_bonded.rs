use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::{ArgonRules, NbList};
use bioshell_ff::Energy;
use bioshell_ff::{Coordinates, System};
use bioshell_numerical::Vec3;
use bioshell_sim::generators::random_chain;

#[test]
fn simple_harmonic_test() {
    // ---------- Create a simple system of 3 atoms
    let mut xyz = Coordinates::new(3);
    for i in 0..3 {
        xyz[i].x = i as f64 * 10.0;
    }
    // ---------- Create system's list of neighbors
    let mut nbl: NbList = NbList::new(2.0, 1.0, Box::new(ArgonRules {}));
    // ---------- Create a system
    let mut system: System = System::new(xyz, nbl);

    // ---------- Harmonic energy
    let springs = SimpleHarmonic::new(10.0, 1.0);
    assert!(f64::abs(springs.energy(&system)) < 0.00001);

    // ---------- Check energy change while moving the last atom away
    let energies = vec![0.25, 1.0, 2.25, 4.0];
    for dx in 0..4 {
        system.add(0, 0.5, 0.0, 0.0);
        assert!(f64::abs(springs.energy(&system) - energies[dx]) < 0.0001);
    }

    // ---------- Create a system of 101 atoms, i.e. 100 springs
    xyz = Coordinates::new(101);
    random_chain(9.5, 1.0, &Vec3::new(0.0, 0.0, 0.0), &mut xyz);

    // ---------- Create system's list of neighbors
    nbl = NbList::new(2.0, 1.0, Box::new(ArgonRules {}));
    system = System::new(xyz, nbl);

    assert!(f64::abs(springs.energy(&system) - 25.0) < 0.01);
}
