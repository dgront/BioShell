use bioshell_ff::nonbonded::{NbList, PairwiseNonbondedEvaluator, PolymerRules, SimpleContact};
use bioshell_ff::Energy;
use bioshell_ff::{Coordinates, System};

#[test]
fn simple_contact_test() {
    // ---------- Create a simple system of 3 atoms
    let mut xyz = Coordinates::new(3);
    for i in 0..3 {
        xyz[i].x = 2.0 * i as f64;
    }
    // ---------- Contact energy  2.3, 3.3,4.3,10.0,-1.0
    let contacts = PairwiseNonbondedEvaluator::new(
        4.5,
        Box::new(SimpleContact::new(2.3, 3.3, 4.3, 10.0, -1.0)),
    );

    // ---------- Create system's list of neighbors
    let nbl: NbList = NbList::new(4.5, 1.0, Box::new(PolymerRules {}));

    // ---------- Create a system, update NB list
    let mut system: System = System::new(xyz, nbl);

    // ---------- For this system all the 3 atoms are within repulsion range, no exclusion rules
    assert!(f64::abs(contacts.energy(&system) + 1.0) < 0.0001);

    // ---------- Check energy change while moving the last atom away
    system.set(2, 2.0, 0.0, 0.0);
    let energies = vec![10.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, 0.0]; // 2.25, 2.5, 2.75, 3, 3.25 ... 4.25, 4.5
    for dx in 0..energies.len() {
        system.add(2, 0.25, 0.0, 0.0);
        // println!("{} {} {}",system.coordinates()[2].x, contacts.energy(&system), energies[dx]);
        assert!(f64::abs(contacts.energy(&system) - energies[dx]) < 0.00001);
    }
}
