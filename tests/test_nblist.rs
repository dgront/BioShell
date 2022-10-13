use rand::Rng;

use bioshell_ff::Coordinates;
use bioshell_ff::nonbonded::{NbList, ArgonRules, PolymerRules};
use bioshell_numerical::Vec3;
use bioshell_sim::generators::cubic_grid_atoms;

macro_rules! check_nb_count {
    ($nbl:expr, $n_atoms:expr, $expected:expr) => {
        for i in 0..$n_atoms {
            assert_eq!($nbl.neighbors(i).len(), $expected);
        }
    }
}

macro_rules! assert_vec3_eq {
    ($vi:expr, $vj:expr, $epsilon:expr) => {
        assert!(f64::abs($vi.x - $vj.x) < $epsilon);
        assert!(f64::abs($vi.y - $vj.y) < $epsilon);
        assert!(f64::abs($vi.z - $vj.z) < $epsilon);
    }
}

#[allow(non_snake_case)]
fn create_1D_system(N:usize) -> Coordinates {
    let mut xyz = Coordinates::new(N);
    xyz.set_box_len(N as f64);
    for i in 0..N {
        xyz[i].x = i as f64 + 0.5;
    }
    xyz
}

#[test]
fn simple_nblist_full_update() {
    // ---------- Create a simple system of N = 5 atoms
    const N : usize = 5;
    let xyz = create_1D_system(N);

    // ---------- Create a list of neighbors for that system
    let mut nbl = NbList::new(0.8,0.5, Box::new(ArgonRules{}));
    assert_eq!(nbl.cutoff(), 0.8);              // --- check cutoff and border size
    assert_eq!(nbl.buffer_width(), 0.5);

    nbl.update_all(&xyz);               // --- update the list
    check_nb_count!(nbl, N, 2);                 // --- check if the number of neighbors is OK

    // --- modify cutoff, buffer and check neighbors again
    nbl.set_cutoff(1.4);
    assert_eq!(nbl.cutoff(), 1.4);
    nbl.update_all(&xyz);
    check_nb_count!(nbl, N, 2);

    nbl.set_buffer_width(0.9);
    assert_eq!(nbl.buffer_width(), 0.9);
    nbl.update_all(&xyz);
    check_nb_count!(nbl, N, 4);
}

#[test]
fn simple_nblist_incremental_update() {
    // ---------- Create a simple system of N = 5 atoms
    const N: usize = 5;
    let mut xyz = create_1D_system(N);
    // ---------- Create a list of neighbors for that system
    let mut nbl = NbList::new(0.8,0.5, Box::new(ArgonRules{}));
    nbl.update_all(&xyz);               // --- build the list

    for _ in 0..10 {
        xyz.add(2, 0.2, 0.0, 0.0);
        nbl.update(&xyz, 2);
    }
}

macro_rules! check_nbl {
    ($coords:expr, $nbl:expr, $N:expr, $d_max:expr) => {
        for i in 0..$N {
            for j in 0..i {
                if $coords.closest_distance_square(i,j).sqrt() <= $d_max {
                    assert!($nbl.neighbors(i).contains(&j));
                    assert!($nbl.neighbors(j).contains(&i));
                }
            }
        }
    }
}


#[test]
fn test_polymer_rules() {

    // ---------- Create a simple system of 3 atoms
    let mut xyz = Coordinates::new(3);
    for i in 0..3 {
        xyz[i].x = 4.0 * i as f64;
    }
    // ---------- Create system's list of neighbors
    let mut nbl: NbList = NbList::new(4.5,1.0,Box::new(PolymerRules{}));
    nbl.update_all(&xyz);

    for i in 0..3 {
        assert_eq!(nbl.neighbors(i).len(), 0);
    }
}

#[test]
fn nblist_3d_test() {

    let max_step: f64 = 0.2;

    // ---------- Create a simple system of N = 5 atoms
    const N: usize = 27;
    let mut coords = Coordinates::new(N);

    // ---------- Distribute atoms on a 3x3x3 grid
    coords.set_box_len(3.0);
    cubic_grid_atoms(&mut coords);
    assert_vec3_eq!(coords[0], Vec3::new(0.5, 0.5, 0.5), 0.000001);
    assert_vec3_eq!(coords[13], Vec3::new(1.5, 1.5, 1.5), 0.000001);

    // ---------- Create a list of neighbors for that system
    let mut nbl = NbList::new(0.6,0.5, Box::new(ArgonRules{}));
    nbl.update_all(&coords);               // --- build the list

    check_nbl!(coords, nbl, N, nbl.cutoff());

    let mut rng = rand::thread_rng();

    for _ in 0..10000 {
        let i_moved = rng.gen_range(0..N);
        coords.add(i_moved,rng.gen_range(-max_step..max_step),
                   rng.gen_range(-max_step..max_step),rng.gen_range(-max_step..max_step));
        nbl.update(&coords, i_moved);

        check_nbl!(coords, nbl, N, nbl.cutoff());
    }
}
