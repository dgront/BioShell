use bioshell_ff::Coordinates;
use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_ff::Energy;
use bioshell_numerical::Vec3;

#[test]
fn simple_contact_test() {
    // ---------- Create a simple system of 3 atoms
    let mut xyz = Coordinates::new(3);
    for i in 0..3 {
        xyz[i].x = i as f32;
    }
    // ---------- Contact energy  2.3, 3.3,4.3,10.0,-1.0
    let contacts = PairwiseNonbondedEvaluator::new(1, 4.5 as f32,
            Box::new(SimpleContact::new(2.3,3.3,4.3,10.0,-1.0)));

    // ---------- Check energy change while moving the last atom away
    let energies = vec![10.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.0, -1.0, -1.0, 0.0]; // 2.25, 2.5, 2.75, 3, 3.25 ... 4.25, 4.5
    for dx in 0..energies.len() {
        xyz[2].x += 0.25;
        // println!("{} {} {}",xyz[2].x, contacts.energy(&xyz), energies[dx]);
        assert!(f64::abs(contacts.energy(&xyz) - energies[dx]) < 0.00001);
    }
}
