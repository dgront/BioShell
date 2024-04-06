use bioshell_pdb::calc::Vec3;
use bioshell_pdb::{load_pdb_file, load_pdb_reader};

use surpass::{HBond3CA, SurpassAlphaSystem, SurpassEnergy};


#[test]
fn test_energy_6ca() {
    let vip = Vec3::new(-3.0, 0.0, 0.0);
    let vi = Vec3::new(0.0, 2.3, 0.0);
    let vin = Vec3::new(3.0, 0.0, 0.0);
    let mut vjp = Vec3::new(-3.0, 0.0, 5.0);
    let mut vj = Vec3::new(0.0, 2.3, 5.0);
    let mut vjn = Vec3::new(3.0, 0.0, 5.0);
    let hb_energy = HBond3CA::new();
    let en = hb_energy.evaluate_hbond_energy_6ca(1, &vip, &vi, &vin, 10, &vjp, &vj, &vjn);
    println!("{}", en);
    assert!(en > 0.99);
    let v_one = Vec3::new(0.0, 0.0, 1.0);
    vjp += &v_one;
    vj += &v_one;
    vjn += &v_one;
    let en = hb_energy.evaluate_hbond_energy_6ca(1, &vip, &vi, &vin, 10, &vjp, &vj, &vjn);
    println!("{}", en);
    assert!(en > 0.99);
    vjn += &v_one;
    let en = hb_energy.evaluate_hbond_energy_6ca(1, &vip, &vi, &vin, 10, &vjp, &vj, &vjn);
    assert!(en < 0.01);
}

#[test]
fn test_energy_2gb1() {
    let model = load_2gb1();
    let hb_energy = HBond3CA::new();
    let en = hb_energy.evaluate_hbond_energy(&model, 43, 52);
    println!("{}", en);
}

#[test]
fn test_surpass_energy() {
    let model = load_2gb1();
    let hb_energy = HBond3CA::new();
    let en = hb_energy.evaluate(&model);
    println!("{}", en);
}

fn load_2gb1() -> SurpassAlphaSystem {
    let str_2gb1 = load_pdb_file("tests/test_inputs/2gb1.pdb").unwrap();
    return  SurpassAlphaSystem::from_pdb_structure(&str_2gb1, 100.0);
}