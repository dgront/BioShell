#[derive(Clone, Debug)]
/// Counts how many system perturbations were successful.
///
/// Each Monte Carlo [`Mover`](Mover) must contain an [AcceptanceStatistics]
/// and update its counters accordingly to the outcome of a [`Mover::perturb()`] call.
/// The total number of Monte Carlo moves attempted is `n_succ + n_failed`
pub struct AcceptanceStatistics {
    /// number of successful perturbations
    pub n_succ: i32,
    /// number of failures
    pub n_failed: i32,
}

impl AcceptanceStatistics {
    /// Computes the success rate for a given Monte Carlo Markov chain.
    ///
    /// Simply returns `n_succ / (n_succ + n_failed)`
    pub fn success_rate(&self) -> f64 {
        let sum = self.n_succ + self.n_failed;
        if sum == 0 {
            return 0.0;
        }
        return self.n_succ as f64 / (sum as f64);
    }

    /// Computes the success rate since the given point in simulation
    ///
    /// The success rate is computed based on *new* observations that were made
    /// after the given `prev_stats` were recorded
    pub fn recent_success_rate(&self, prev_stats: &AcceptanceStatistics) -> f64 {
        let succ = self.n_succ - prev_stats.n_succ;
        let fail = self.n_failed - prev_stats.n_failed;
        let sum = succ + fail;
        if sum == 0 {
            return 0.0;
        }
        return succ as f64 / (sum as f64);
    }
}

impl Default for AcceptanceStatistics {
    fn default() -> Self {
        AcceptanceStatistics {
            n_succ: 0,
            n_failed: 0,
        }
    }
}