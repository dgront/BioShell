use rand::Rng;
use std::ops::Range;
use crate::CartesianSystem;
use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Energy, System};
use bioshell_numerical::vec3::Vec3;
use bioshell_numerical::{random_point_nearby, Rototranslation};


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
