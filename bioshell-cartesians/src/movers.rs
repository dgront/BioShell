use rand::Rng;
use std::ops::Range;
use crate::CartesianSystem;
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Energy, System};
use bioshell_numerical::vec3::Vec3;
use bioshell_numerical::{random_point_nearby, Rototranslation};

/// A mover that moves a single, randomly selected atom by a small random vector.
pub struct SingleAtomMove {
    max_step: f64,
    succ_rate: AcceptanceStatistics,
}

impl SingleAtomMove {
    /// Create a new mover that shifts a single atom.
    pub fn new(max_range: f64) -> SingleAtomMove {
        SingleAtomMove {
            max_step: max_range,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for SingleAtomMove {
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
        let mut rng = rand::thread_rng();//obtain a random number generator
        let i_moved = rng.gen_range(0..system.size());//obtain a random number

        let old_x: f64 = system.coordinates()[i_moved].x;//obtain a random coordinate of x
        let old_y: f64 = system.coordinates()[i_moved].y;//obtain a random coordinate of y
        let old_z: f64 = system.coordinates()[i_moved].z;//obtain a random coordinate of z
        let old_energy: f64 = energy.energy_by_pos(system, i_moved);

        //add the randomly generated coordinate into the Cartesian System.
        system.add(
            i_moved,
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
            rng.gen_range(-self.max_step..self.max_step),
        );

        let new_energy: f64 = energy.energy_by_pos(system, i_moved);

        if acc.check(old_energy, new_energy)
        {
            //success
            system.update_nbl(i_moved);//update non-bonded list of neighbors
            self.succ_rate.n_succ += 1;
            return Option::from(i_moved..i_moved);
        }
        else
        {
            //failure
            self.succ_rate.n_failed += 1;
            system.set(i_moved, old_x, old_y, old_z);//fall back
            return Option::None;
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_step
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_step = new_val;
    }
}

/// A mover that changes a volume of a Cartesian system.
pub struct ChangeVolume {
    /// pressure of the system `$p$`
    pub pressure: f64,
    pub temperature: f64,
    max_step: f64,
    succ_rate: AcceptanceStatistics,
}

#[allow(non_upper_case_globals)]
const pV_to_Kelvins: f64 = 1.0E-30 / 1.380649E-23; // 10^-30 divided by the Boltzmann constant

impl ChangeVolume {
    pub fn new(pressure: f64, temperature: f64) -> ChangeVolume {
        ChangeVolume {
            pressure,
            temperature,
            max_step: 0.01,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for ChangeVolume {
    #[allow(non_snake_case)]
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        _acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
        // ---------- attempt the volume change
        let en_before = energy.energy(system);
        let mut rng = rand::thread_rng();
        let v0 = system.volume();
        let lnV0 = v0.ln();
        let lnV = lnV0 + rng.gen_range(-self.max_step..self.max_step);
        let new_V = lnV.exp();
        let old_len = system.box_len();
        let new_len = new_V.powf(0.333333333);
        let f = new_len / system.box_len();
        for i in 0..system.size() {
            let x = system.coordinates().x(i) * f;
            let y = system.coordinates().y(i) * f;
            let z = system.coordinates().z(i) * f;
            system.set(i, x, y, z);
        }
        system.set_box_len(new_len);
        let en_after = energy.energy(system);
        let weight = ((en_after - en_before) + self.pressure * (new_V - v0) * pV_to_Kelvins)
            / self.temperature
            - (system.size() + 1) as f64 * (new_V / v0).ln();

        if rng.gen_range(0.0..1.0) > (-weight / self.temperature).exp() {
            // --- move rejected
            for i in 0..system.size() {
                let x = system.coordinates().x(i) / f;
                let y = system.coordinates().y(i) / f;
                let z = system.coordinates().z(i) / f;
                system.set(i, x, y, z);
            }
            self.succ_rate.n_failed += 1;
            system.set_box_len(old_len);
            return None;
        } else {
            self.succ_rate.n_succ += 1;
            return Some(0..system.size());
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_step
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_step = new_val;
    }
}

/// A mover that applies a crankshaft move to a molecule.
pub struct CrankshaftMove {
    max_angle: f64,
    frag_size: usize,//how far would be the axis positions
    succ_rate: AcceptanceStatistics,
}

impl CrankshaftMove {
    /// Create a new mover that applies a crankshaft move.
    pub fn new(max_angle: f64) -> CrankshaftMove {
        CrankshaftMove {
            max_angle,
            frag_size:5,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for CrankshaftMove {
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>>
    {
        let mut rng = rand::thread_rng();//obtain a random number generator
        let system_length = system.size();//obtain the length of the Cartesian system
        let i_start = rng.gen_range(0..system_length - &self.frag_size-1);//obtain a random start position
        let mut i_end = i_start + &self.frag_size + 1;//obtain the end position relative to the start position
        let angle = rng.gen_range(-&self.max_angle..self.max_angle);//obtain a random angle
        let mut start = system.coordinates()[i_start].clone();
        let mut end = system.coordinates()[i_end].clone();
        let roto_tran = Rototranslation::around_axis(&start, &end, angle);

        let energy_before = energy.energy(system);//calculate the energy before the move

        //apply forward rotation
        for i in i_start +1..i_end  {
            let mut temp_coord:Vec3 = system.coordinates()[i].clone();
            //print!("Before ({},{},{})", temp_coord.x, temp_coord.y, temp_coord.z);
            roto_tran.apply_mut(&mut temp_coord);
            //println!("After ({},{},{})", temp_coord.x, temp_coord.y, temp_coord.z);
            system.set_vec(i, temp_coord);
        }

        let energy_after = energy.energy(system);

        if acc.check(energy_before, energy_after)
        {
            //if move succeeds
            self.succ_rate.n_succ += 1;
            return Option::from(i_start..i_end + 1);
        }
        else
        {
            //if move doesn't succeed
            self.succ_rate.n_failed += 1;
            //fall back - apply inverse rotation
            for i in i_start +1..i_end  {
                let mut temp_coord:Vec3 = system.coordinates()[i].clone();
                roto_tran.apply_inverse_mut(&mut temp_coord);
                system.set_vec(i, temp_coord);
            }
            return None;
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_angle
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_angle = new_val;
    }
}



/// A mover that applies a terminal move to a molecule.
pub struct TerminalMove {
    max_angle: f64,
    frag_size: usize,
    succ_rate: AcceptanceStatistics,
}

// The move involves a rotation of one or two peptide bonds at one of the chain termini picked
// randomly with rotation axis passing through the alpha carbon in a random direction.
// It is important that in this scheme any rotation and its inverse are picked with equal
// probability.
impl TerminalMove {
    /// Create a new mover that applies a terminal move.
    pub fn new(max_angle: f64) -> TerminalMove {
        TerminalMove {
            max_angle,
            frag_size:5,
            succ_rate: Default::default(),
        }
    }
}

impl<E: Energy<CartesianSystem>> Mover<CartesianSystem, E> for TerminalMove {
    fn perturb(
        &mut self,
        system: &mut CartesianSystem,
        energy: &E,
        acc: &mut dyn AcceptanceCriterion,
    ) -> Option<Range<usize>> {
        let mut rng = rand::thread_rng();
        let size = system.size();
        let terminal_flag = rng.gen_range(0..1 + 1);

        let mut i_start = 0;
        let mut i_end = 0;
        let mut terminal_vec = Vec3::zero();
        let mut term_rot_axis_start = Vec3::zero();

        if terminal_flag == 0 {
            i_start = 0;
            i_end = self.frag_size;
            terminal_vec = system.coordinates()[i_end-1];
            term_rot_axis_start = system.coordinates()[i_end];
        } else {
            i_start = size - self.frag_size;
            i_end = size;
            term_rot_axis_start = system.coordinates()[i_start];
            terminal_vec = system.coordinates()[i_start+1];
        }

        let angle = rng.gen_range(-&self.max_angle..self.max_angle);
        let start = term_rot_axis_start;
        // let end = random_point_nearby(&start, self.frag_size as f64);
        let end = terminal_vec;
        let roto_tran = Rototranslation::around_axis(&start, &end, angle);
        let energy_before = energy.energy(system);//calculate the energy before the move

        // --- apply forward rotation
        for i in i_start..i_end {
            let mut temp_coord: Vec3 = system.coordinates()[i].clone();
            roto_tran.apply_mut(&mut temp_coord);
            system.set_vec(i, temp_coord);
        }

        let energy_after = energy.energy(system);

        if acc.check(energy_before, energy_after) { // --- if move succeeds
            self.succ_rate.n_succ += 1;
            return Option::from(i_start..i_end + 1);
        } else {            // --- if move doesn't succeed
            self.succ_rate.n_failed += 1;
            // --- fall back - apply inverse rotation
            for i in i_start..i_end {
                let mut temp_coord: Vec3 = system.coordinates()[i].clone();
                roto_tran.apply_inverse_mut(&mut temp_coord);
                system.set_vec(i, temp_coord);
            }
            return None;
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics {
        self.succ_rate.clone()
    }

    fn max_range(&self) -> f64 {
        self.max_angle
    }

    fn set_max_range(&mut self, new_val: f64) {
        self.max_angle = new_val;
    }
}
