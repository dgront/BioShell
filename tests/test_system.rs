use bioshell_ff::{Coordinates, System};
use bioshell_ff::nonbonded::{ArgonRules, NbList};
use bioshell_sim::generators::cubic_grid_atoms;

#[test]
fn create_simple_system() {
    const L: f32 = 6.0;        // --- width of the box
    const N: usize = 27;       // --- the number of atoms in it
    const R: f32 = 2.0;         // --- radius of each atom (interaction cutoff)
    const B: f32 = 1.0;         // --- buffer zone width for NBL

    // ---------- Create system's coordinates and initialize them
    let mut xyz = Coordinates::new(N);
    xyz.set_box_len(L);
    cubic_grid_atoms(&mut xyz);

    // ---------- Create system's list of neighbors - it's mandatory!
    let nbl: NbList = NbList::new(2.0*R,B,Box::new(ArgonRules{}));

    // ---------- Create and modify the system
    let mut system: System = System::new(xyz,nbl);
    system.add(0, 0.5, 0.5, 0.5);

}