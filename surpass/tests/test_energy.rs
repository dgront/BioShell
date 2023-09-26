use surpass::{SurpassAlphaSystem, ExcludedVolume, SurpassEnergy, MoveProposal};

macro_rules! assert_delta {
    ($x:expr, $y:expr, $d:expr) => {
        assert!(($x-$y).abs() < $d, "a = {}, b = {}", $x, $y)
    }
}


/// Test ExcludedVolume energy for a system of 5 atoms where atom 2 clashes with atoms 0 and 4
#[test]
fn test_excluded_volume_5() {
    // --- make a system of 5 atoms
    let mut model = SurpassAlphaSystem::new(&[5], 100.0);
    // --- put all atoms in 0.0,0
    for i in 0..5 {
        model.cax[i] = model.real_to_int(0.0);
        model.cay[i] = model.real_to_int(0.0);
        model.caz[i] = model.real_to_int(0.0);
    }
    // --- move atoms 1 and 3 far away so they don't interact
    model.cax[1] = model.real_to_int(40.0);
    model.cax[3] = model.real_to_int(-40.0);

    // --- set the stage with atoms 4 and 2 - just beyond the interaction cutoff
    model.cax[4] = model.real_to_int(6.0);
    model.cax[2] = model.real_to_int(3.0);
    model.cay[2] = model.real_to_int(4.01);

    // --- propose a move that causes two clashes
    let mut mp: MoveProposal<1> = MoveProposal::new();
    mp.first_moved_pos = 2;
    mp.cay[0] = model.real_to_int(3.99);
    mp.cax[0] = model.real_to_int(3.00);

    let energy = ExcludedVolume::new(&model, 5.0, 100.0);
    assert_delta!(energy.evaluate_delta::<1>(&model, &mp), 200.0, 0.00001);

    assert_delta!(energy.evaluate(&model), 0.0, 0.00001);
    // --- move the atom 2 to make two clashes
    model.cay[2] = model.real_to_int(3.99);
    assert_delta!(energy.evaluate(&model), 200.0, 0.00001);
}

/// Test ExcludedVolume energy where every atom clashes with any other; this should result in (N-2)x(N-1)/2 clashes
#[test]
#[allow(non_snake_case)]
fn test_excluded_volume_N() {

    run_test_excluded_volume_N::<4>();
    run_test_excluded_volume_N::<10>();
    run_test_excluded_volume_N::<13>();
}

#[allow(non_snake_case)]
fn run_test_excluded_volume_N<const N: usize>() {
    // --- make a system of 5 atoms
    let mut model = SurpassAlphaSystem::new(&[N], 100.0);
    // --- put all atoms in 0.0,0
    for i in 0..N {
        model.cax[i] = model.real_to_int(0.0);
        model.cay[i] = model.real_to_int(0.0);
        model.caz[i] = model.real_to_int(0.0);
    }

    let energy = ExcludedVolume::new(&model, 5.0, 100.0);
    assert_delta!(energy.evaluate(&model), ((N-2)*(N-1)) as f64 * 100.0 / 2.0, 0.00001);
}