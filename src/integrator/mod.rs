mod ias15;
mod leapfrog;
pub mod output;
pub mod whfast;

pub use self::ias15::*;
pub use self::leapfrog::*;
pub use self::whfast::WHFast;

use std::any::Any;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

pub trait Integrator {
    fn as_any(&self) -> &dyn Any;
    fn get_n_historic_snapshots(&self) -> usize;
    fn get_n_particles(&self) -> usize;
    fn get_current_time(&self) -> f64;
    fn set_time_limit(&mut self, time_limit: f64);
    fn set_snapshot_periods(
        &mut self,
        historic_snapshot_period: f64,
        recovery_snapshot_period: f64,
    );
    fn initialize_physical_values(&mut self);
    fn iterate(
        &mut self,
        universe_history_writer: &mut BufWriter<File>,
        silent_mode: bool,
    ) -> Result<bool, String>;
    fn write_recovery_snapshot(
        &mut self,
        snapshot_path: &Path,
        universe_history_writer: &mut BufWriter<File>,
    );
}
