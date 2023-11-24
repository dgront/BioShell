use rand::rngs::SmallRng;
use rand::SeedableRng;
use surpass::{SurpassAlphaSystem, ExcludedVolume, SurpassEnergy, MoveProposal, CaContactEnergy, NonBondedEnergy, HingeMove, Mover, NonBondedEnergyKernel};

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
        model.bbx[i] = model.real_to_int(0.0);
        model.bby[i] = model.real_to_int(0.0);
        model.bbz[i] = model.real_to_int(0.0);
    }
    // --- move atoms 1 and 3 far away so they don't interact
    model.bbx[1] = model.real_to_int(40.0);
    model.bbx[3] = model.real_to_int(-40.0);

    // --- set the stage with atoms 4 and 2 - just beyond the interaction cutoff
    model.bbx[4] = model.real_to_int(6.0);
    model.bbx[2] = model.real_to_int(3.0);
    model.bby[2] = model.real_to_int(4.01);

    // --- propose a move that causes two clashes
    let mut mp: MoveProposal<1> = MoveProposal::new();
    mp.first_moved_pos = 2;
    mp.cay[0] = model.real_to_int(3.99);
    mp.cax[0] = model.real_to_int(3.00);

    let excl_vol_kernel = ExcludedVolume::new(&model, 5.0, 100.0);
    let energy: NonBondedEnergy<ExcludedVolume> = NonBondedEnergy::new(&model, excl_vol_kernel);

    assert_delta!(energy.evaluate_delta::<1>(&model, &mp), 200.0, 0.00001);

    assert_delta!(energy.evaluate(&model), 0.0, 0.00001);
    // --- move the atom 2 to make two clashes
    model.bby[2] = model.real_to_int(3.99);
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
    // --- make a system of N atoms
    let mut model = SurpassAlphaSystem::new(&[N], 100.0);
    // --- put all atoms in 0, 0, 0
    for i in 0..N {
        model.bbx[i] = model.real_to_int(0.0);
        model.bby[i] = model.real_to_int(0.0);
        model.bbz[i] = model.real_to_int(0.0);
    }

    let excl_vol_kernel = ExcludedVolume::new(&model, 5.0, 100.0);
    let energy: NonBondedEnergy<ExcludedVolume> = NonBondedEnergy::new(&model,  excl_vol_kernel);
    assert_delta!(energy.evaluate(&model), ((N)*(N-1)) as f64 * 100.0 / 2.0, 0.00001);
}

#[test]
fn test_cacontacts_kernel() {
    let mut system = SurpassAlphaSystem::new(&[10], 1000.0);
    let cntcts = CaContactEnergy::new(&system, 10.0, -1.0, 3.7, 4.0, 5.0);
    let d = [3.2, 3.5, 3.6, 3.701, 3.99, 4.001, 4.999, 5.001];
    let e = [10.0, 10.0, 10.0, 0.0, 0.0, -1.0, -1.0, 0.0];
    for i in 0..d.len() {
        let di = system.real_to_int(d[i]) as f64;
        let en = cntcts.energy_for_distance_squared(di*di);
        assert_delta!(en, e[i], 0.000001);
    }
}

#[test]
fn check_contact_delta_en() {
    const N: usize = 10;
    let mut rnd = SmallRng::seed_from_u64(42);
    let mut system = SurpassAlphaSystem::make_random(&[N], 1000.0, &mut rnd);
    let cntcts = CaContactEnergy::new(&system, 10.0, -1.0, 3.7, 4.0, 5.0);
    let energy: NonBondedEnergy<CaContactEnergy> = NonBondedEnergy::new(&system, cntcts);
    let hinge_mover: HingeMove<4> = HingeMove::new(std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    let mut hinge_prop: MoveProposal<4> = MoveProposal::new();

    let mut en_before: f64; // --- for debugging
    let mut backup : MoveProposal<4> = MoveProposal::new();
    for i in 0..10000 {
        hinge_mover.propose(&mut system, &mut rnd, &mut hinge_prop);
        backup.first_moved_pos = hinge_prop.first_moved_pos;
        backup.backup(&system);
        // eprintln!("--------------------- BEFORE ------------");
        en_before = energy.evaluate(&system);
        let delta_e = energy.evaluate_delta(&system, &hinge_prop);
        hinge_prop.apply(&mut system);
        // eprintln!("--------------------- AFTER ------------");
        let en_after = energy.evaluate(&system);
        if (en_after-en_before - delta_e).abs() > 0.00001 {
            panic!("Incorrect energy evaluation after {i} moves");
        }
    }
}