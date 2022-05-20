use bioshell_ff::Coordinates;
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::Energy;
use bioshell_numerical::Vec3;
use bioshell_sim::generators::random_chain;

#[test]
fn simple_harmonic_test() {
    // ---------- Create a simple system of 3 atoms
    let mut xyz = Coordinates::new(3);
    for i in 0..3 {
        xyz[i].x = i as f32 * 10.0;
    }
    // ---------- Harmonic energy
    let springs = SimpleHarmonic::new(10.0, 1.0);
    assert!(f64::abs(springs.energy(&xyz))<0.00001);

    // ---------- Check energy change while moving the last atom away
    let energies = vec![1.0, 4.0, 9.0, 16.0];
    for dx in 0..4 {
        xyz[2].x += 1.0;
        assert!(f64::abs(springs.energy(&xyz) - energies[dx]) < 0.00001);
    }

    // ---------- Create a system of 101 atoms, i.e. 100 springs
    xyz = Coordinates::new(101);
    random_chain(9.5,1.0,&mut xyz);
    assert!(f64::abs(springs.energy(&xyz)-25.0)<0.00001);
}
