use super::super::Particle;
use super::super::constants::{
    DBL_EPSILON_2, G, IMPLICIT_MIDPOINT_MAX_ITER, IMPLICIT_MIDPOINT_MIN_ITER, MAX_PARTICLES, PI,
    WHFAST_NMAX_NEWT, WHFAST_NMAX_QUART,
};
use super::super::effects::EvolutionType;
use super::super::effects::GeneralRelativityImplementation;
use super::super::particles::Axes;
use super::super::particles::IgnoreGravityTerms;
use super::super::particles::Universe;
use super::Integrator;
use super::output::{write_historic_snapshot, write_recovery_snapshot};
use serde::{Deserialize, Serialize};
use serde_big_array::BigArray;
use std;
use std::any::Any;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::iter;
use std::path::Path;
use time;
use time::{OffsetDateTime, format_description};

/// Source: Rein & Tamayo, 2015
/// `WHFast` is a complete reimplementation of the Wisdom-Holman integrator[1],
/// designed to speed up the algorithm and increase its accuracy (e.g., related
/// to finite double floating-point precision, all real numbers cannot be
/// represented exactly and floating-point precision has consequences for the
/// numerical stability of any algorithm and the growth of numerical round-off
/// error). It is unbiased (i.e., the errors are random and uncorrelated), and
/// has a very slow error growth. For sufficiently small timesteps, it achieves
/// Brouwer’s law (i.e., the energy error grows as time to the power of one
/// half).
///
/// `WHFast` is a symplectic integrator[2], althought its symplectic nature is
/// formally lost as soon as self-gravity or collisions are approximated or
/// when velocity dependent forces are included (such as tidal forces).
///
/// This implementation is equivalent to the REBOUND code (Rein & Liu, 2011) in
/// safe mode (required by the tidal effects) and without correction (i.e. 2nd
/// order integrator, comparable to mercury symplectic part of the hybrid
/// integrator)
///
/// Possible coordinates:
///
/// - Jacobi coordinates: it leads to a better precision if orbits are well
/// separated and do not cross each other.
/// - Democratic-Heliocentric coordinates (a.k.a. canonical heliocentric,
/// Poincare or mixed-variables coordinates): better precision if there are
/// close orbits and/or encounters.
/// - WHDS: it splits the Hamiltonian into three parts (Hernandez and Dehnen,
/// 2017) and solves the two body problem exactly, contrary to the democratic
/// heliocentric one. It keeps the democratic heliocentric properties of being
/// more accurate with close orbits and encounters.
///
/// Sources:
/// - Hernandez and Dehnen, 2017
///     <http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1612.05329>
/// - Ari Silburt, Hanno Rein & Dan Tamayo
///     HERMES: a hybrid integrator for simulating close encounters and planetesimal migration
///     <https://silburt.github.io/files/HERMES.pdf>
///
///
/// [1] Source: Rein & Liu, 2011
///
/// Symplectic Wisdom-Holman mappings (such as the implemented in the SWIFT
/// code4) are mixed variable integrator that calculates the Keplerian motion
/// of two bodies orbiting each other exactly up to machine precision during the
/// drift sub-step. Thus, it is very accurate for problems in which
/// the particle motion is dominated by a central 1/r potential and
/// perturbations added in the kick sub-step are small. Kepler’s equation is
/// solved iteratively every time-step for every particle.
///
///
/// [2] Source: <https://en.wikipedia.org/wiki/N-body_problem#Few_bodies>
///
/// For N > 2, the N-body problem is chaotic,[37] which means that even small errors in
/// integration may grow exponentially in time. Third, a simulation may be over large stretches of
/// model time (e.g. millions of years) and numerical errors accumulate as integration time
/// increases.
///
/// There are a number of techniques to reduce errors in numerical integration.[18] Local
/// coordinate systems are used to deal with widely differing scales in some problems, for example
/// an Earth-Moon coordinate system in the context of a solar system simulation. Variational
/// methods and perturbation theory can yield approximate analytic trajectories upon which the
/// numerical integration can be a correction. The use of a symplectic integrator ensures that the
/// simulation obeys Hamilton's equations to a high degree of accuracy and in particular that
/// energy is conserved.
///

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum CoordinatesType {
    Jacobi,
    DemocraticHeliocentric,
    WHDS,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct WHFast {
    time_step: f64,
    half_time_step: f64,
    pub universe: Universe,
    pub current_time: f64,
    current_iteration: usize,
    pub recovery_snapshot_period: f64,
    pub historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
    /// Internal data structures below. Nothing to be changed by the user.
    #[serde(with = "BigArray")]
    particles_alternative_coordinates: [AlternativeCoordinates; MAX_PARTICLES], // Jacobi, democractic-heliocentric or WHDS
    alternative_coordinates_type: CoordinatesType,
    timestep_warning: usize,
    #[serde(with = "BigArray")]
    inertial_velocity_errors: [Axes; MAX_PARTICLES], // A running compensation for lost low-order bits (Kahan 1965; Higham 2002; Hairer et al. 2006)
    #[serde(with = "BigArray")]
    particle_angular_momentum_errors: [Axes; MAX_PARTICLES], // A running compensation for lost low-order bits (Kahan 1965; Higham 2002; Hairer et al. 2006)
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct AlternativeCoordinates {
    mass: f64,
    mass_g: f64,
    // Jacobi, democractic-heliocentric or WHDS
    position: Axes,
    velocity: Axes,
    acceleration: Axes,
}

impl Hash for WHFast {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Hash only works with certain types (for instance, it does not work with f64 values)
        // thus we convert the whole integrator to a string thanks to the debug trait
        // and we hash that value
        format!("{self:?}").hash(state);
    }
}

impl WHFast {
    pub fn new(
        time_step: f64,
        recovery_snapshot_period: f64,
        historic_snapshot_period: f64,
        universe: Universe,
        alternative_coordinates_type: CoordinatesType,
    ) -> WHFast {
        let particles_alternative_coordinates = [AlternativeCoordinates {
            mass: 0.,
            mass_g: 0.,
            position: Axes::new(),
            velocity: Axes::new(),
            acceleration: Axes::new(),
        }; MAX_PARTICLES];

        WHFast {
            time_step,
            half_time_step: 0.5 * time_step,
            recovery_snapshot_period,
            historic_snapshot_period,
            last_recovery_snapshot_time: -1.,
            last_historic_snapshot_time: -1.,
            n_historic_snapshots: 0,
            hash: 0,
            universe,
            current_time: 0.,
            current_iteration: 0,
            // WHFast specifics:
            particles_alternative_coordinates,
            alternative_coordinates_type,
            timestep_warning: 0,
            inertial_velocity_errors: [Axes::new(); MAX_PARTICLES],
            particle_angular_momentum_errors: [Axes::new(); MAX_PARTICLES],
        }
    }
}

impl Integrator for WHFast {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn get_n_historic_snapshots(&self) -> usize {
        self.n_historic_snapshots
    }

    fn get_n_particles(&self) -> usize {
        self.universe.n_particles
    }

    fn get_current_time(&self) -> f64 {
        self.current_time
    }

    fn set_time_limit(&mut self, time_limit: f64) {
        if time_limit > 0. && self.universe.time_limit != time_limit {
            if time_limit > self.universe.time_limit && self.universe.consider_effects.evolution {
                // Check if the new time is in the range of the evolutionary model
                for (i, evolver) in self.universe.particles_evolvers.iter().enumerate() {
                    // If it is an evolving body.
                    if !matches!(evolver.evolution, EvolutionType::NonEvolving)
                        && evolver.time[evolver.time.len() - 1] < time_limit
                    {
                        panic!(
                            "Your new time limit ({} days) is greater than the maximum allowed age of the evolving body #{} ({} days)",
                            time_limit,
                            i + 1,
                            evolver.time[evolver.time.len() - 1]
                        );
                    }
                }
            } else if time_limit < self.universe.time_limit && time_limit < self.current_time {
                panic!(
                    "Your new time limit ({} days) is smaller than the current time ({} days)",
                    time_limit, self.current_time
                );
            }
            println!(
                "[INFO {} UTC] The time limit changed from {} to {} days",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                self.universe.time_limit,
                time_limit
            );
            self.universe.time_limit = time_limit;
        }
    }

    fn set_snapshot_periods(
        &mut self,
        historic_snapshot_period: f64,
        recovery_snapshot_period: f64,
    ) {
        if historic_snapshot_period > 0.
            && self.historic_snapshot_period != historic_snapshot_period
        {
            println!(
                "[INFO {} UTC] The historic snapshot period changed from {} to {} days",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                self.historic_snapshot_period,
                historic_snapshot_period
            );
            self.historic_snapshot_period = historic_snapshot_period;
        } else {
            println!(
                "[INFO {} UTC] A historic snapshot will be saved every {} days",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                self.historic_snapshot_period
            );
        }

        if recovery_snapshot_period > 0.
            && self.recovery_snapshot_period != recovery_snapshot_period
        {
            println!(
                "[INFO {} UTC] The recovery snapshot period changed from {} to {} days",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                self.recovery_snapshot_period,
                recovery_snapshot_period
            );
            self.recovery_snapshot_period = recovery_snapshot_period;
        } else {
            println!(
                "[INFO {} UTC] A recovery snapshot will be saved every {} days",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                self.recovery_snapshot_period
            );
        }
    }

    fn initialize_physical_values(&mut self) {
        assert!(
            !(self.current_time != 0.),
            "Physical values cannot be initialized on a resumed simulation"
        );
        let evolution = true;
        self.universe
            .calculate_spin_and_evolving_quantities(self.current_time, evolution); // Make sure we start with the good initial values
        self.universe.calculate_roche_radiuses(); // Needed for collision detection
    }

    fn iterate(
        &mut self,
        universe_history_writer: &mut BufWriter<File>,
        silent_mode: bool,
    ) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger =
            self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger =
            self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            self.universe.inertial_to_heliocentric();
            let evolution = true;
            self.universe
                .calculate_spin_and_evolving_quantities(self.current_time, evolution);
            if self.universe.consider_effects.tides {
                self.universe.calculate_denergy_dt();
            }
            write_historic_snapshot(
                universe_history_writer,
                &self.universe,
                self.current_time,
                self.time_step,
            );
            if first_snapshot_trigger {
                self.last_historic_snapshot_time = 0.;
            } else {
                // Do not use `self.current_time` to avoid small deviations
                // Do not use `self.n_historic_snapshots as f64*self.historic_snapshot_period` because `historic_snapshot_period` can be changed by the user when resuming an already started simulation
                self.last_historic_snapshot_time += self.historic_snapshot_period;
            }
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time / 365.25;
            if !silent_mode {
                print!(
                    "Year: {:0.0} ({:0.1e}) | Time step: {:0.3} days                    \r",
                    current_time_years, current_time_years, self.time_step
                );
                let _ = std::io::stdout().flush();
            }
        }

        let ignored_gravity_terms = match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => IgnoreGravityTerms::WHFastOne,
            CoordinatesType::WHDS | CoordinatesType::DemocraticHeliocentric => {
                IgnoreGravityTerms::WHFastTwo
            }
        };
        let ignore_gravity_terms = ignored_gravity_terms;
        //
        let mut consider_dangular_momentum_dt_from_general_relativity = false;
        if self.universe.consider_effects.general_relativity
            && self.universe.general_relativity_implementation
                == GeneralRelativityImplementation::Kidder1995
        {
            consider_dangular_momentum_dt_from_general_relativity = true;
        }
        let integrate_spin = self.universe.consider_effects.tides
            || self.universe.consider_effects.rotational_flattening
            || self.universe.consider_effects.evolution
            || consider_dangular_momentum_dt_from_general_relativity;

        let evolution = true;
        self.integrate_velocity_dependent_forces(self.half_time_step, integrate_spin, evolution); // Corrects the inertial velocity and computes spin
        // A 'DKD'-like integrator will do the first 'D' part:
        self.iterate_position_and_velocity_with_whfasthelio_part1(); // updates alternative pos/vel (not using the inertial acceleration but the keplerian motion) by half-step
        self.universe
            .gravity_calculate_acceleration(ignore_gravity_terms); // computes gravitational inertial accelerations
        //
        //// Posidonius' additional effects are velocity dependent, thus we do not compute them here:
        //self.universe.inertial_to_heliocentric(); // required to compute additional effects
        //let evolution = true;
        //let dangular_momentum_dt = false;
        //let accelerations = true;
        //self.universe.calculate_additional_effects(self.current_time, evolution, dangular_momentum_dt, accelerations, ignored_gravity_terms); // changes inertial accelerations
        //
        // A 'DKD'-like integrator will do the 'KD' part:
        self.iterate_position_and_velocity_with_whfasthelio_part2(); // updates alternative and inertial pos/vel (using the inertial acceleration by one time step to compute interactions and the keplerian motion by half-step)
        let evolution = false; // Only evolve once per full step (optimization)
        self.integrate_velocity_dependent_forces(self.half_time_step, integrate_spin, evolution); // Corrects the inertial velocity and computes spin
        self.current_time += self.time_step;

        // ---------------------------------------------------------------------
        self.current_iteration += 1;

        // Return
        if self.current_time + self.time_step > self.universe.time_limit {
            Err("Simulation completed".to_string())
        } else {
            Ok(first_snapshot_trigger || recovery_snapshot_time_trigger)
        }
    }

    fn write_recovery_snapshot(
        &mut self,
        snapshot_path: &Path,
        universe_history_writer: &mut BufWriter<File>,
    ) {
        self.last_recovery_snapshot_time = self.current_time;
        universe_history_writer.flush().unwrap();
        // Compute hash for this universe at this moment of time
        self.hash = 0;
        let mut s = DefaultHasher::new();
        self.hash(&mut s);
        self.hash = s.finish();
        write_recovery_snapshot(snapshot_path, &self);
    }
}

impl WHFast {
    fn integrate_velocity_dependent_forces(
        &mut self,
        dt: f64,
        integrate_spin: bool,
        evolution: bool,
    ) {
        // Integrate velocity-dependent forces using an implicit midpoint
        // - It follows the schema implemented in REBOUNDx and described in Tamayo et al. 2020
        // - It corresponds to the implicit midpoint method (https://en.wikipedia.org/wiki/Midpoint_method)
        // - In addition, it implements compensated summation
        let ignored_gravity_terms = match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => IgnoreGravityTerms::WHFastOne,
            CoordinatesType::WHDS | CoordinatesType::DemocraticHeliocentric => {
                IgnoreGravityTerms::WHFastTwo
            }
        };

        let mut particles_orig = self.universe.particles;
        let mut particles_final: [Particle; MAX_PARTICLES] = particles_orig;
        let mut particles_prev: [Particle; MAX_PARTICLES];
        let mut inertial_velocity_changes: [Axes; MAX_PARTICLES] = [Axes::new(); MAX_PARTICLES];
        let mut angular_momentum_changes: [Axes; MAX_PARTICLES] = [Axes::new(); MAX_PARTICLES];
        let mut converged = false;
        for i in 0..IMPLICIT_MIDPOINT_MAX_ITER {
            particles_prev = particles_final;
            // To calculate non-gravity/additional accelerations:
            // - Positions and velocities are needed in heliocentric
            // - But additional accelerations are computed in inertial (i.e., barycentric)
            self.universe.inertial_to_heliocentric(); // required to compute additional effects
            // Calculate non-gravity accelerations
            {
                let evolution = evolution && i == 0; // Only evolve in the first iteration (optimization)
                let dangular_momentum_dt = integrate_spin;
                let accelerations = true;
                self.universe.calculate_additional_effects(
                    self.current_time,
                    evolution,
                    dangular_momentum_dt,
                    accelerations,
                    ignored_gravity_terms,
                ); // changes inertial accelerations
            }
            // Compute final velocity/spin angular momentum
            for (
                (
                    (
                        (((particle_avg, particle_orig), particle_final), inertial_velocity_error),
                        angular_momentum_error,
                    ),
                    inertial_velocity_change,
                ),
                angular_momentum_change,
            ) in self.universe.particles[..self.universe.n_particles]
                .iter()
                .zip(particles_orig[..self.universe.n_particles].iter_mut())
                .zip(particles_final[..self.universe.n_particles].iter_mut())
                .zip(self.inertial_velocity_errors[..self.universe.n_particles].iter())
                .zip(self.particle_angular_momentum_errors[..self.universe.n_particles].iter())
                .zip(inertial_velocity_changes[..self.universe.n_particles].iter_mut())
                .zip(angular_momentum_changes[..self.universe.n_particles].iter_mut())
            {
                // Compensated summation to improve the accuracy of additions that involve
                // one small and one large floating point number (already considered by IAS15)
                //      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
                //      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
                inertial_velocity_change.x = dt * particle_avg.inertial_additional_acceleration.x
                    - inertial_velocity_error.x;
                inertial_velocity_change.y = dt * particle_avg.inertial_additional_acceleration.y
                    - inertial_velocity_error.y;
                inertial_velocity_change.z = dt * particle_avg.inertial_additional_acceleration.z
                    - inertial_velocity_error.z;
                //
                particle_final.inertial_velocity.x =
                    particle_orig.inertial_velocity.x + inertial_velocity_change.x;
                particle_final.inertial_velocity.y =
                    particle_orig.inertial_velocity.y + inertial_velocity_change.y;
                particle_final.inertial_velocity.z =
                    particle_orig.inertial_velocity.z + inertial_velocity_change.z;
                if integrate_spin {
                    // Compensated summation to improve the accuracy of additions that involve
                    // one small and one large floating point number (already considered by IAS15)
                    //      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
                    //      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
                    angular_momentum_change.x =
                        dt * particle_avg.dangular_momentum_dt.x - angular_momentum_error.x;
                    angular_momentum_change.y =
                        dt * particle_avg.dangular_momentum_dt.y - angular_momentum_error.y;
                    angular_momentum_change.z =
                        dt * particle_avg.dangular_momentum_dt.z - angular_momentum_error.z;
                    //
                    particle_final.angular_momentum.x =
                        particle_orig.angular_momentum.x + angular_momentum_change.x;
                    particle_final.angular_momentum.y =
                        particle_orig.angular_momentum.y + angular_momentum_change.y;
                    particle_final.angular_momentum.z =
                        particle_orig.angular_momentum.z + angular_momentum_change.z;
                }
            }
            // Compare final with previous velocity/spin but make sure there is a minimum of
            // iterations first to guarantee a minimum precision independent of the time step
            if i >= IMPLICIT_MIDPOINT_MIN_ITER - 1
                && self.converged_velocity_dependent_forces_integration(
                    &particles_final,
                    &particles_prev,
                    integrate_spin,
                )
            {
                converged = true;
                break;
            }
            // Average velocities and spins using original and final ones
            // - Updates velocity/spin in self.universe.particles (i.e., particles_avg)
            self.average_particles_for_velocity_dependent_forces_integration(
                &particles_orig,
                &particles_final,
                integrate_spin,
            );
        }
        if !converged {
            println!(
                "[WARNING {} UTC] WHFast convergence issue with the integration of the additional forces. Most probably the perturbation is too strong.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        //
        for (
            (
                (
                    (((particle, particle_final), particle_orig), inertial_velocity_error),
                    angular_momentum_error,
                ),
                inertial_velocity_change,
            ),
            angular_momentum_change,
        ) in self.universe.particles[..self.universe.n_particles]
            .iter_mut()
            .zip(particles_final[..self.universe.n_particles].iter())
            .zip(particles_orig[..self.universe.n_particles].iter())
            .zip(self.inertial_velocity_errors[..self.universe.n_particles].iter_mut())
            .zip(self.particle_angular_momentum_errors[..self.universe.n_particles].iter_mut())
            .zip(inertial_velocity_changes[..self.universe.n_particles].iter())
            .zip(angular_momentum_changes[..self.universe.n_particles].iter_mut())
        {
            // Inertial positions and accelerations didn't change
            // Use modified velocities
            particle.inertial_velocity.x = particle_final.inertial_velocity.x;
            particle.inertial_velocity.y = particle_final.inertial_velocity.y;
            particle.inertial_velocity.z = particle_final.inertial_velocity.z;

            // Compensated summation to improve the accuracy of additions that involve
            // one small and one large floating point number (already considered by IAS15)
            //      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
            //      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
            inertial_velocity_error.x = (particle.inertial_velocity.x
                - particle_orig.inertial_velocity.x)
                - inertial_velocity_change.x;
            inertial_velocity_error.y = (particle.inertial_velocity.y
                - particle_orig.inertial_velocity.y)
                - inertial_velocity_change.y;
            inertial_velocity_error.z = (particle.inertial_velocity.z
                - particle_orig.inertial_velocity.z)
                - inertial_velocity_change.z;
            if integrate_spin {
                // Use modified angular_momentum
                particle.angular_momentum.x = particle_final.angular_momentum.x;
                particle.angular_momentum.y = particle_final.angular_momentum.y;
                particle.angular_momentum.z = particle_final.angular_momentum.z;

                // Compensated summation to improve the accuracy of additions that involve
                // one small and one large floating point number (already considered by IAS15)
                //      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
                //      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
                angular_momentum_error.x = (particle.angular_momentum.x
                    - particle_orig.angular_momentum.x)
                    - angular_momentum_change.x;
                angular_momentum_error.y = (particle.angular_momentum.y
                    - particle_orig.angular_momentum.y)
                    - angular_momentum_change.y;
                angular_momentum_error.z = (particle.angular_momentum.z
                    - particle_orig.angular_momentum.z)
                    - angular_momentum_change.z;
            }
        }
    }

    fn converged_velocity_dependent_forces_integration(
        &mut self,
        particles_final: &[Particle; MAX_PARTICLES],
        particles_prev: &[Particle; MAX_PARTICLES],
        integrate_spin: bool,
    ) -> bool {
        let mut final_total_velocity_2 = 0.;
        let mut delta_total_velocity_2 = 0.;
        let mut final_total_angular_momentum_2 = 0.;
        let mut delta_total_angular_momentum_2 = 0.;
        for (particle_final, particle_prev) in particles_final[..self.universe.n_particles]
            .iter()
            .zip(particles_prev[..self.universe.n_particles].iter())
        {
            // Velocity
            let dvx = particle_final.inertial_velocity.x - particle_prev.inertial_velocity.x;
            let dvy = particle_final.inertial_velocity.y - particle_prev.inertial_velocity.y;
            let dvz = particle_final.inertial_velocity.z - particle_prev.inertial_velocity.z;
            delta_total_velocity_2 += dvx.powi(2) + dvy.powi(2) + dvz.powi(2);
            final_total_velocity_2 += particle_final.inertial_velocity.x.powi(2)
                + particle_final.inertial_velocity.y.powi(2)
                + particle_final.inertial_velocity.z.powi(2);
            if integrate_spin {
                // Angular momentum
                let dsx = particle_final.angular_momentum.x - particle_prev.angular_momentum.x;
                let dsy = particle_final.angular_momentum.y - particle_prev.angular_momentum.y;
                let dsz = particle_final.angular_momentum.z - particle_prev.angular_momentum.z;
                delta_total_angular_momentum_2 += dsx.powi(2) + dsy.powi(2) + dsz.powi(2);
                final_total_angular_momentum_2 += particle_final.angular_momentum.x.powi(2)
                    + particle_final.angular_momentum.y.powi(2)
                    + particle_final.angular_momentum.z.powi(2);
            }
        }

        (delta_total_angular_momentum_2 / final_total_angular_momentum_2 < DBL_EPSILON_2
            || !integrate_spin)
            && delta_total_velocity_2 / final_total_velocity_2 < DBL_EPSILON_2
    }

    fn average_particles_for_velocity_dependent_forces_integration(
        &mut self,
        particles_orig: &[Particle; MAX_PARTICLES],
        particles_final: &[Particle; MAX_PARTICLES],
        integrate_spin: bool,
    ) {
        for ((particle_avg, particle_orig), particle_final) in self.universe.particles
            [..self.universe.n_particles]
            .iter_mut()
            .zip(particles_orig[..self.universe.n_particles].iter())
            .zip(particles_final[..self.universe.n_particles].iter())
        {
            // Velocity
            particle_avg.inertial_velocity.x =
                0.5 * (particle_orig.inertial_velocity.x + particle_final.inertial_velocity.x);
            particle_avg.inertial_velocity.y =
                0.5 * (particle_orig.inertial_velocity.y + particle_final.inertial_velocity.y);
            particle_avg.inertial_velocity.z =
                0.5 * (particle_orig.inertial_velocity.z + particle_final.inertial_velocity.z);
            if integrate_spin {
                // Angular momentum
                particle_avg.angular_momentum.x =
                    0.5 * (particle_orig.angular_momentum.x + particle_final.angular_momentum.x);
                particle_avg.angular_momentum.y =
                    0.5 * (particle_orig.angular_momentum.y + particle_final.angular_momentum.y);
                particle_avg.angular_momentum.z =
                    0.5 * (particle_orig.angular_momentum.z + particle_final.angular_momentum.z);
            }
        }
    }

    // WHFast integrator
    fn iterate_position_and_velocity_with_whfasthelio_part1(&mut self) {
        let half_time_step = self.half_time_step;

        // ---------------------------------------------------------------------
        // A 'DKD'-like integrator will do the first 'D' part.
        self.inertial_to_alternative_posvel();
        self.kepler_steps(half_time_step);
        self.jump_step(half_time_step);
        self.alternative_to_inertial_posvel();
        // ---------------------------------------------------------------------
    }

    fn iterate_position_and_velocity_with_whfasthelio_part2(&mut self) {
        let time_step = self.time_step;
        let half_time_step = self.half_time_step;

        // Continue the WHFAST integration ('KD' part)
        self.interaction_step(time_step); // changes alternative velocities using inertial accelerations
        self.jump_step(half_time_step);
        self.kepler_steps(half_time_step);
        self.alternative_to_inertial_posvel();
    }

    //------------------------------
    //Operators
    fn jump_step(&mut self, dt: f64) {
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => Self::jacobi_jump_step(dt),
            CoordinatesType::DemocraticHeliocentric => self.democratic_heliocentric_jump_step(dt),
            CoordinatesType::WHDS => self.whds_jump_step(dt),
        }
    }

    fn jacobi_jump_step(_dt: f64) {
        // Nothing to be done
    }

    fn democratic_heliocentric_jump_step(&mut self, dt: f64) {
        let (particles_left, particles_right) = self.universe.particles
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, particles_right)) = particles_right.split_first_mut() {
            let m0 = star.mass;
            let mut px = 0.;
            let mut py = 0.;
            let mut pz = 0.;
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((_, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    px += particle.mass * particle_alternative_coordinates.velocity.x;
                    py += particle.mass * particle_alternative_coordinates.velocity.y;
                    pz += particle.mass * particle_alternative_coordinates.velocity.z;
                }
                for particle_alternative_coordinates in particles_alternative_coordinates_left
                    .iter_mut()
                    .chain(particles_alternative_coordinates_right.iter_mut())
                {
                    particle_alternative_coordinates.position.x += dt * px / m0;
                    particle_alternative_coordinates.position.y += dt * py / m0;
                    particle_alternative_coordinates.position.z += dt * pz / m0;
                }
            }
        }
    }

    fn whds_jump_step(&mut self, dt: f64) {
        let (particles_left, particles_right) = self.universe.particles
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, particles_right)) = particles_right.split_first_mut() {
            let m0 = star.mass;
            let mut px = 0.;
            let mut py = 0.;
            let mut pz = 0.;
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((_, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    let factor = m0 + particle.mass; // mf numerator
                    px += particle.mass * particle_alternative_coordinates.velocity.x / factor;
                    py += particle.mass * particle_alternative_coordinates.velocity.y / factor;
                    pz += particle.mass * particle_alternative_coordinates.velocity.z / factor;
                }
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    let factor = m0 + particle.mass; // mf numerator
                    particle_alternative_coordinates.position.x += dt
                        * (px
                            - (particle.mass * particle_alternative_coordinates.velocity.x
                                / factor));
                    particle_alternative_coordinates.position.y += dt
                        * (py
                            - (particle.mass * particle_alternative_coordinates.velocity.y
                                / factor));
                    particle_alternative_coordinates.position.z += dt
                        * (pz
                            - (particle.mass * particle_alternative_coordinates.velocity.z
                                / factor));
                }
            }
        }
    }

    fn interaction_step(&mut self, dt: f64) {
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_interaction_step(dt),
            CoordinatesType::DemocraticHeliocentric => {
                self.democratic_heliocentric_interaction_step(dt);
            }
            CoordinatesType::WHDS => self.whds_interaction_step(dt),
        }
    }

    fn jacobi_interaction_step(&mut self, dt: f64) {
        self.inertial_to_jacobi_acc();
        let softening = 1e-12;
        let (_, particles_right) = self.universe.particles[..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, _)) = particles_right.split_first_mut() {
            let m0 = star.mass;
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((_, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                let mut eta = m0;
                for (i, particle_alternative_coordinates) in particles_alternative_coordinates_left
                    .iter_mut()
                    .chain(particles_alternative_coordinates_right.iter_mut())
                    .enumerate()
                {
                    eta += particle_alternative_coordinates.mass;
                    particle_alternative_coordinates.velocity.x +=
                        dt * particle_alternative_coordinates.acceleration.x;
                    particle_alternative_coordinates.velocity.y +=
                        dt * particle_alternative_coordinates.acceleration.y;
                    particle_alternative_coordinates.velocity.z +=
                        dt * particle_alternative_coordinates.acceleration.z;
                    if i > 0 {
                        // The star is already out of the loop, ignore first planet too (it should be the closer to the star)
                        let rj2i = 1.
                            / (particle_alternative_coordinates.position.x.powi(2)
                                + particle_alternative_coordinates.position.y.powi(2)
                                + particle_alternative_coordinates.position.z.powi(2)
                                + softening);
                        let rji = sqrt!(rj2i);
                        let rj3im = rji * rj2i * G * eta;
                        let prefac = dt * rj3im;
                        particle_alternative_coordinates.velocity.x +=
                            prefac * particle_alternative_coordinates.position.x;
                        particle_alternative_coordinates.velocity.y +=
                            prefac * particle_alternative_coordinates.position.y;
                        particle_alternative_coordinates.velocity.z +=
                            prefac * particle_alternative_coordinates.position.z;
                    }
                }
            }
        }
    }

    fn democratic_heliocentric_interaction_step(&mut self, dt: f64) {
        let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
            self.particles_alternative_coordinates[..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((_, particles_alternative_coordinates_right)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            let (particles_left, particles_right) = self.universe.particles
                [..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((_, particles_right)) = particles_right.split_first_mut() {
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    particle_alternative_coordinates.velocity.x +=
                        dt * particle.inertial_acceleration.x;
                    particle_alternative_coordinates.velocity.y +=
                        dt * particle.inertial_acceleration.y;
                    particle_alternative_coordinates.velocity.z +=
                        dt * particle.inertial_acceleration.z;
                }
            }
        }
    }

    fn whds_interaction_step(&mut self, dt: f64) {
        let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
            self.particles_alternative_coordinates[..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((_, particles_alternative_coordinates_right)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            let (particles_left, particles_right) = self.universe.particles
                [..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star, particles_right)) = particles_right.split_first_mut() {
                let m0 = star.mass;
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    let factor = m0 + particle.mass; // mf numerator
                    particle_alternative_coordinates.velocity.x +=
                        dt * factor * particle.inertial_acceleration.x / m0;
                    particle_alternative_coordinates.velocity.y +=
                        dt * factor * particle.inertial_acceleration.y / m0;
                    particle_alternative_coordinates.velocity.z +=
                        dt * factor * particle.inertial_acceleration.z / m0;
                }
            }
        }
    }

    fn kepler_steps(&mut self, time_step: f64) {
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_kepler_steps(time_step),
            CoordinatesType::DemocraticHeliocentric => {
                self.democratic_heliocentric_kepler_steps(time_step);
            }
            CoordinatesType::WHDS => self.whds_kepler_steps(time_step),
        }
        let (_, particles_alternative_coordinates_right) = self.particles_alternative_coordinates
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star_alternative_coordinates, _)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            star_alternative_coordinates.position.x +=
                time_step * star_alternative_coordinates.velocity.x;
            star_alternative_coordinates.position.y +=
                time_step * star_alternative_coordinates.velocity.y;
            star_alternative_coordinates.position.z +=
                time_step * star_alternative_coordinates.velocity.z;
        }
    }

    fn jacobi_kepler_steps(&mut self, time_step: f64) {
        let mut star_planets_mass_g =
            self.universe.particles[self.universe.hosts.index.most_massive].mass_g;
        for i in 0..self.universe.n_particles {
            if i == self.universe.hosts.index.most_massive {
                continue;
            }
            star_planets_mass_g += self.particles_alternative_coordinates[i].mass_g;
            self.kepler_individual_step(i, star_planets_mass_g, time_step);
        }
    }

    fn democratic_heliocentric_kepler_steps(&mut self, time_step: f64) {
        let star_mass_g = self.universe.particles[self.universe.hosts.index.most_massive].mass_g;
        for i in 0..self.universe.n_particles {
            if i == self.universe.hosts.index.most_massive {
                continue;
            }
            self.kepler_individual_step(i, star_mass_g, time_step);
        }
    }

    fn whds_kepler_steps(&mut self, time_step: f64) {
        let star_mass_g = self.universe.particles[self.universe.hosts.index.most_massive].mass_g;
        for i in 0..self.universe.n_particles {
            if i == self.universe.hosts.index.most_massive {
                continue;
            }
            let star_planet_mass_g = star_mass_g + self.particles_alternative_coordinates[i].mass_g;
            self.kepler_individual_step(i, star_planet_mass_g, time_step);
        }
    }

    //-----------------------------
    // Keplerian motion
    fn kepler_individual_step(&mut self, i: usize, mass_g: f64, dt: f64) {
        // Save a copy of the original position and velocities
        let p1_position;
        let p1_velocity;
        {
            let p1 = &self.particles_alternative_coordinates[i];
            p1_position = p1.position;
            p1_velocity = p1.velocity;
        }

        let r0 = sqrt!(p1_position.x.powi(2) + p1_position.y.powi(2) + p1_position.z.powi(2));
        let r0i = 1. / r0;
        let v2 = p1_velocity.x.powi(2) + p1_velocity.y.powi(2) + p1_velocity.z.powi(2);
        let beta = 2. * mass_g * r0i - v2;
        let eta0 = p1_position.x * p1_velocity.x
            + p1_position.y * p1_velocity.y
            + p1_position.z * p1_velocity.z;
        let zeta0 = mass_g - beta * r0;
        let mut x;
        let mut gs;
        let mut invperiod = 0.; // only used for beta>0.
        let x_per_period;

        if beta > 0. {
            //// Elliptic orbit
            let sqrt_beta = sqrt!(beta);
            invperiod = sqrt_beta * beta / (2. * PI * mass_g);
            x_per_period = 2. * PI / sqrt_beta;
            if abs!(dt) * invperiod > 1. && self.timestep_warning == 0 {
                // Ignoring const qualifiers. This warning should not have any effect on
                // other parts of the code, nor is it vital to show it.
                println!(
                    "[WARNING {} UTC] WHFast convergence issue. Timestep is larger than at least one orbital period.",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap()
                );
                self.timestep_warning += 1;
            }
            //x =  dt*invperiod*x_per_period; // first order guess
            let dtr0i = dt * r0i;
            //x = dtr0i; // first order guess
            x = dtr0i * (1. - dtr0i * eta0 * 0.5 * r0i); // second order guess
        //x = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
        //x =  dt*beta/mass_g + eta0/mass_g*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version
        } else {
            //// Hyperbolic orbit
            x = 0.; // Initial guess
            x_per_period = f64::NAN; // only used for beta>0. nan triggers Newton's method for beta<0
        }

        let mut converged = 0;
        let mut old_x = x;

        //// Do one Newton step
        gs = WHFast::stiefel_gs3(beta, x);
        let eta0_gs1_zeta0_gs2 = eta0 * gs[1] + zeta0 * gs[2];
        let mut ri = 1. / (r0 + eta0_gs1_zeta0_gs2);
        x = ri * (x * eta0_gs1_zeta0_gs2 - eta0 * gs[2] - zeta0 * gs[3] + dt);

        // Choose solver depending on estimated step size
        // Note, for hyperbolic orbits this uses Newton's method.
        if !x_per_period.is_nan() && abs!(x - old_x) > 0.01 * x_per_period {
            // Quartic solver
            // Linear initial guess
            x = beta * dt / mass_g;
            let mut prev_x = [0.; WHFAST_NMAX_QUART + 1];
            for n_lag in 1..WHFAST_NMAX_QUART {
                gs = WHFast::stiefel_gs3(beta, x);
                let f = r0 * x + eta0 * gs[2] + zeta0 * gs[3] - dt;
                let fp = r0 + eta0 * gs[1] + zeta0 * gs[2];
                let fpp = eta0 * gs[0] + zeta0 * gs[1];
                let denom = fp + sqrt!(abs!(16. * fp * fp - 20. * f * fpp));
                x = (x * denom - 5. * f) / denom;
                if prev_x.iter().any(|&x_prev_i| x == x_prev_i) {
                    // Converged. Exit.
                    //n_lag = WHFAST_NMAX_QUART;
                    converged = 1;
                    break;
                }
                prev_x[n_lag] = x;
            }
            let eta0_gs1_zeta0_gs2 = eta0 * gs[1] + zeta0 * gs[2];
            ri = 1. / (r0 + eta0_gs1_zeta0_gs2);
        } else {
            // Newton's method
            let mut old_x2;
            for _ in 1..WHFAST_NMAX_NEWT {
                old_x2 = old_x;
                old_x = x;
                gs = WHFast::stiefel_gs3(beta, x);
                let eta0_gs1_zeta0_gs2 = eta0 * gs[1] + zeta0 * gs[2];
                ri = 1. / (r0 + eta0_gs1_zeta0_gs2);
                x = ri * (x * eta0_gs1_zeta0_gs2 - eta0 * gs[2] - zeta0 * gs[3] + dt);

                if x == old_x || x == old_x2 {
                    // Converged. Exit.
                    converged = 1;
                    break;
                }
            }
        }

        // If solver did not work, fallback to bisection
        if converged == 0 {
            let mut x_min;
            let mut x_max;
            if beta > 0. {
                //Elliptic
                x_min = x_per_period * (dt * invperiod).floor();
                x_max = x_min + x_per_period;
            } else {
                //Hyperbolic
                let h2 = r0 * r0 * v2 - eta0 * eta0;
                let q = h2 / mass_g / (1. + sqrt!(1. - h2 * beta / (mass_g * mass_g)));
                let vq = sqrt!(h2) / q;
                x_min = 1. / (vq + r0 / dt);
                x_max = dt / q;
            }
            x = (x_max + x_min) / 2.;
            loop {
                gs = WHFast::stiefel_gs3(beta, x);
                let s = r0 * x + eta0 * gs[2] + zeta0 * gs[3] - dt;
                if s >= 0. {
                    x_max = x;
                } else {
                    x_min = x;
                }
                x = (x_max + x_min) / 2.;

                if abs!(x_max - x_min) / x_max <= 1e-15 {
                    break;
                }
            }
            let eta0_gs1_zeta0_gs2 = eta0 * gs[1] + zeta0 * gs[2];
            ri = 1. / (r0 + eta0_gs1_zeta0_gs2);
        }
        if ri.is_nan() {
            // Exception for (almost) straight line motion in hyperbolic case
            ri = 0.;
            gs[1] = 0.;
            gs[2] = 0.;
            gs[3] = 0.;
        }

        // Note: These are not the traditional f and g functions.
        let f = -mass_g * gs[2] * r0i;
        let g = dt - mass_g * gs[3];
        let fd = -mass_g * gs[1] * r0i * ri;
        let gd = -mass_g * gs[2] * ri;

        let p_j = &mut self.particles_alternative_coordinates[i];
        p_j.position.x += f * p1_position.x + g * p1_velocity.x;
        p_j.position.y += f * p1_position.y + g * p1_velocity.y;
        p_j.position.z += f * p1_position.z + g * p1_velocity.z;

        // WARNING: p1_position used below should not be modified by the previous block (keep p_j
        // independent)
        p_j.velocity.x += fd * p1_position.x + gd * p1_velocity.x;
        p_j.velocity.y += fd * p1_position.y + gd * p1_velocity.y;
        p_j.velocity.z += fd * p1_position.z + gd * p1_velocity.z;
    }

    fn stiefel_gs3(beta: f64, x: f64) -> [f64; 6] {
        let x2 = x.powi(2);
        let mut gs = WHFast::stumpff_cs3(beta * x2);
        gs[1] *= x;
        gs[2] *= x2;
        gs[3] *= x2 * x;
        gs
    }

    fn stumpff_cs3(z: f64) -> [f64; 6] {
        // Fast inverse factorial lookup table
        let invfactorial: [f64; 35] = [
            1.,
            1.,
            1. / 2.,
            1. / 6.,
            1. / 24.,
            1. / 120.,
            1. / 720.,
            1. / 5040.,
            1. / 40_320.,
            1. / 362_880.,
            1. / 3_628_800.,
            1. / 39_916_800.,
            1. / 479_001_600.,
            1. / 6_227_020_800.,
            1. / 87_178_291_200.,
            1. / 1_307_674_368_000.,
            1. / 20_922_789_888_000.,
            1. / 355_687_428_096_000.,
            1. / 6_402_373_705_728_000.,
            1. / 121_645_100_408_832_000.,
            1. / 2_432_902_008_176_640_000.,
            1. / 51_090_942_171_709_440_000.,
            1. / 1_124_000_727_777_607_680_000.,
            1. / 25_852_016_738_884_976_640_000.,
            1. / 620_448_401_733_239_439_360_000.,
            1. / 15_511_210_043_330_985_984_000_000.,
            1. / 403_291_461_126_605_635_584_000_000.,
            1. / 10_888_869_450_418_352_160_768_000_000.,
            1. / 304_888_344_611_713_860_501_504_000_000.,
            1. / 8_841_761_993_739_701_954_543_616_000_000.,
            1. / 265_252_859_812_191_058_636_308_480_000_000.,
            1. / 8_222_838_654_177_922_817_725_562_880_000_000.,
            1. / 263_130_836_933_693_530_167_218_012_160_000_000.,
            1. / 8_683_317_618_811_886_495_518_194_401_280_000_000.,
            1. / 295_232_799_039_604_140_847_618_609_643_520_000_000.,
        ];
        let mut z = z;
        let mut n = 0;
        while abs!(z) > 0.1 {
            z /= 4.;
            n += 1;
        }
        let mut cs = [0.; 6];
        let nmax = 13;
        let mut c_odd = invfactorial[nmax];
        let mut c_even = invfactorial[nmax - 1];

        let mut np = nmax - 2;
        while np >= 3 {
            c_odd = invfactorial[np] - z * c_odd;
            c_even = invfactorial[np - 1] - z * c_even;
            np -= 2;
        }
        cs[3] = c_odd;
        cs[2] = c_even;
        cs[1] = invfactorial[1] - z * c_odd;
        cs[0] = invfactorial[0] - z * c_even;
        while n > 0 {
            cs[3] = (cs[2] + cs[0] * cs[3]) * 0.25;
            cs[2] = cs[1] * cs[1] * 0.5;
            cs[1] *= cs[0];
            cs[0] = 2. * cs[0] * cs[0] - 1.;
            n -= 1;
        }
        cs
    }

    //-----------------------------
    // Coordinate transformations
    //-----------------------------
    fn inertial_to_alternative_posvel(&mut self) {
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.inertial_to_jacobi_posvel(),
            CoordinatesType::DemocraticHeliocentric => {
                self.inertial_to_democratic_heliocentric_posvel();
            }
            CoordinatesType::WHDS => self.inertial_to_whds_posvel(),
        }
    }

    fn inertial_to_jacobi_posvel(&mut self) {
        let (particles_left, particles_right) = self.universe.particles
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, particles_right)) = particles_right.split_first_mut() {
            let m0 = star.mass;
            let mut eta = m0;
            let mut s_x = eta * star.inertial_position.x;
            let mut s_y = eta * star.inertial_position.y;
            let mut s_z = eta * star.inertial_position.z;
            let mut s_vx = eta * star.inertial_velocity.x;
            let mut s_vy = eta * star.inertial_velocity.y;
            let mut s_vz = eta * star.inertial_velocity.z;
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    let ei = 1. / eta;
                    eta += particle.mass;
                    let pme = eta * ei;
                    particle_alternative_coordinates.mass = particle.mass;
                    particle_alternative_coordinates.mass_g = particle.mass_g;
                    particle_alternative_coordinates.position.x =
                        particle.inertial_position.x - s_x * ei;
                    particle_alternative_coordinates.position.y =
                        particle.inertial_position.y - s_y * ei;
                    particle_alternative_coordinates.position.z =
                        particle.inertial_position.z - s_z * ei;
                    particle_alternative_coordinates.velocity.x =
                        particle.inertial_velocity.x - s_vx * ei;
                    particle_alternative_coordinates.velocity.y =
                        particle.inertial_velocity.y - s_vy * ei;
                    particle_alternative_coordinates.velocity.z =
                        particle.inertial_velocity.z - s_vz * ei;
                    s_x = s_x * pme + particle.mass * particle_alternative_coordinates.position.x;
                    s_y = s_y * pme + particle.mass * particle_alternative_coordinates.position.y;
                    s_z = s_z * pme + particle.mass * particle_alternative_coordinates.position.z;
                    s_vx = s_vx * pme + particle.mass * particle_alternative_coordinates.velocity.x;
                    s_vy = s_vy * pme + particle.mass * particle_alternative_coordinates.velocity.y;
                    s_vz = s_vz * pme + particle.mass * particle_alternative_coordinates.velocity.z;
                }
                let mtot = eta;
                let mtot_i = 1. / mtot;
                star_alternative_coordinates.mass = mtot;
                star_alternative_coordinates.position.x = s_x * mtot_i;
                star_alternative_coordinates.position.y = s_y * mtot_i;
                star_alternative_coordinates.position.z = s_z * mtot_i;
                star_alternative_coordinates.velocity.x = s_vx * mtot_i;
                star_alternative_coordinates.velocity.y = s_vy * mtot_i;
                star_alternative_coordinates.velocity.z = s_vz * mtot_i;
            }
        }
    }

    fn inertial_to_jacobi_acc(&mut self) {
        let (particles_left, particles_right) = self.universe.particles
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, particles_right)) = particles_right.split_first_mut() {
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                let mut eta = star.mass;
                let mut s_ax = eta * star.inertial_acceleration.x;
                let mut s_ay = eta * star.inertial_acceleration.y;
                let mut s_az = eta * star.inertial_acceleration.z;
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter().chain(particles_right.iter()))
                {
                    let ei = 1. / eta;
                    eta += particle.mass;
                    let pme = eta * ei;
                    particle_alternative_coordinates.acceleration.x =
                        particle.inertial_acceleration.x - s_ax * ei;
                    particle_alternative_coordinates.acceleration.y =
                        particle.inertial_acceleration.y - s_ay * ei;
                    particle_alternative_coordinates.acceleration.z =
                        particle.inertial_acceleration.z - s_az * ei;
                    s_ax = s_ax * pme
                        + particle.mass * particle_alternative_coordinates.acceleration.x;
                    s_ay = s_ay * pme
                        + particle.mass * particle_alternative_coordinates.acceleration.y;
                    s_az = s_az * pme
                        + particle.mass * particle_alternative_coordinates.acceleration.z;
                }
                let mtot = eta;
                let mtot_i = 1. / mtot;
                star_alternative_coordinates.acceleration.x = s_ax * mtot_i;
                star_alternative_coordinates.acceleration.y = s_ay * mtot_i;
                star_alternative_coordinates.acceleration.z = s_az * mtot_i;
            }
        }
    }

    fn inertial_to_democratic_heliocentric_posvel(&mut self) {
        self.inertial_to_whds_and_democratic_heliocentric_posvel();
    }

    fn inertial_to_whds_posvel(&mut self) {
        self.inertial_to_whds_and_democratic_heliocentric_posvel();
    }

    fn inertial_to_whds_and_democratic_heliocentric_posvel(&mut self) {
        let (particles_left, particles_right) = self.universe.particles
            [..self.universe.n_particles]
            .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star, particles_right)) = particles_right.split_first_mut() {
            let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
                self.particles_alternative_coordinates[..self.universe.n_particles]
                    .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
                particles_alternative_coordinates_right.split_first_mut()
            {
                star_alternative_coordinates.position.zero();
                star_alternative_coordinates.velocity.zero();
                star_alternative_coordinates.mass = 0.;
                star_alternative_coordinates.mass_g = 0.;
                for particle in iter::once(&*star)
                    .chain(particles_left.iter())
                    .chain(particles_right.iter())
                {
                    star_alternative_coordinates.position.x +=
                        particle.inertial_position.x * particle.mass;
                    star_alternative_coordinates.position.y +=
                        particle.inertial_position.y * particle.mass;
                    star_alternative_coordinates.position.z +=
                        particle.inertial_position.z * particle.mass;
                    star_alternative_coordinates.velocity.x +=
                        particle.inertial_velocity.x * particle.mass;
                    star_alternative_coordinates.velocity.y +=
                        particle.inertial_velocity.y * particle.mass;
                    star_alternative_coordinates.velocity.z +=
                        particle.inertial_velocity.z * particle.mass;
                    star_alternative_coordinates.mass += particle.mass;
                    star_alternative_coordinates.mass_g += particle.mass_g;
                }
                let mtot = star_alternative_coordinates.mass;
                star_alternative_coordinates.position.x /= mtot;
                star_alternative_coordinates.position.y /= mtot;
                star_alternative_coordinates.position.z /= mtot;
                star_alternative_coordinates.velocity.x /= mtot;
                star_alternative_coordinates.velocity.y /= mtot;
                star_alternative_coordinates.velocity.z /= mtot;

                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter().chain(particles_right.iter()))
                {
                    particle_alternative_coordinates.position.x =
                        particle.inertial_position.x - star.inertial_position.x;
                    particle_alternative_coordinates.position.y =
                        particle.inertial_position.y - star.inertial_position.y;
                    particle_alternative_coordinates.position.z =
                        particle.inertial_position.z - star.inertial_position.z;
                    particle_alternative_coordinates.velocity.x =
                        particle.inertial_velocity.x - star_alternative_coordinates.velocity.x;
                    particle_alternative_coordinates.velocity.y =
                        particle.inertial_velocity.y - star_alternative_coordinates.velocity.y;
                    particle_alternative_coordinates.velocity.z =
                        particle.inertial_velocity.z - star_alternative_coordinates.velocity.z;
                    particle_alternative_coordinates.mass = particle.mass;
                    particle_alternative_coordinates.mass_g = particle.mass_g;
                    if self.alternative_coordinates_type == CoordinatesType::WHDS {
                        let factor = (star.mass + particle.mass) / star.mass; // mf complete
                        particle_alternative_coordinates.velocity.x *= factor;
                        particle_alternative_coordinates.velocity.y *= factor;
                        particle_alternative_coordinates.velocity.z *= factor;
                    }
                }
            }
        }
    }

    fn alternative_to_inertial_posvel(&mut self) {
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_to_inertial_posvel(),
            CoordinatesType::DemocraticHeliocentric => {
                self.democratic_heliocentric_to_inertial_posvel();
            }
            CoordinatesType::WHDS => self.whds_to_inertial_posvel(),
        }
    }

    fn jacobi_to_inertial_posvel(&mut self) {
        let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
            self.particles_alternative_coordinates[..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            let m0 = star_alternative_coordinates.mass;
            let mut eta = m0;
            let mut s_x = eta * star_alternative_coordinates.position.x;
            let mut s_y = eta * star_alternative_coordinates.position.y;
            let mut s_z = eta * star_alternative_coordinates.position.z;
            let mut s_vx = eta * star_alternative_coordinates.velocity.x;
            let mut s_vy = eta * star_alternative_coordinates.velocity.y;
            let mut s_vz = eta * star_alternative_coordinates.velocity.z;
            let (particles_left, particles_right) = self.universe.particles
                [..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star, particles_right)) = particles_right.split_first_mut() {
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter()
                        .chain(particles_alternative_coordinates_right.iter())
                        .rev()
                        .zip(
                            particles_left
                                .iter_mut()
                                .chain(particles_right.iter_mut())
                                .rev(),
                        )
                {
                    let ei = 1. / eta;
                    s_x = (s_x - particle.mass * particle_alternative_coordinates.position.x) * ei;
                    s_y = (s_y - particle.mass * particle_alternative_coordinates.position.y) * ei;
                    s_z = (s_z - particle.mass * particle_alternative_coordinates.position.z) * ei;
                    s_vx =
                        (s_vx - particle.mass * particle_alternative_coordinates.velocity.x) * ei;
                    s_vy =
                        (s_vy - particle.mass * particle_alternative_coordinates.velocity.y) * ei;
                    s_vz =
                        (s_vz - particle.mass * particle_alternative_coordinates.velocity.z) * ei;
                    particle.inertial_position.x =
                        particle_alternative_coordinates.position.x + s_x;
                    particle.inertial_position.y =
                        particle_alternative_coordinates.position.y + s_y;
                    particle.inertial_position.z =
                        particle_alternative_coordinates.position.z + s_z;
                    particle.inertial_velocity.x =
                        particle_alternative_coordinates.velocity.x + s_vx;
                    particle.inertial_velocity.y =
                        particle_alternative_coordinates.velocity.y + s_vy;
                    particle.inertial_velocity.z =
                        particle_alternative_coordinates.velocity.z + s_vz;
                    eta -= particle.mass;
                    s_x *= eta;
                    s_y *= eta;
                    s_z *= eta;
                    s_vx *= eta;
                    s_vy *= eta;
                    s_vz *= eta;
                }
                let mtot = eta;
                let mtot_i = 1. / mtot;
                star.inertial_position.x = s_x * mtot_i;
                star.inertial_position.y = s_y * mtot_i;
                star.inertial_position.z = s_z * mtot_i;
                star.inertial_velocity.x = s_vx * mtot_i;
                star.inertial_velocity.y = s_vy * mtot_i;
                star.inertial_velocity.z = s_vz * mtot_i;
            }
        }
    }

    fn democratic_heliocentric_to_inertial_posvel(&mut self) {
        self.whds_and_democratic_heliocentric_to_inertial_posvel();
    }

    fn whds_to_inertial_posvel(&mut self) {
        self.whds_and_democratic_heliocentric_to_inertial_posvel();
    }

    fn whds_and_democratic_heliocentric_to_inertial_posvel(&mut self) {
        self.whds_and_democratic_heliocentric_to_inertial_pos();
        let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
            self.particles_alternative_coordinates[..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            let (particles_left, particles_right) = self.universe.particles
                [..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star, particles_right)) = particles_right.split_first_mut() {
                let m0 = star.mass;
                let mut new_star_velocity = star_alternative_coordinates.velocity;
                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    if self.alternative_coordinates_type == CoordinatesType::WHDS {
                        let factor = (m0 + particle.mass) / m0; // mf complete
                        particle.inertial_velocity.x = particle_alternative_coordinates.velocity.x
                            / factor
                            + star_alternative_coordinates.velocity.x;
                        particle.inertial_velocity.y = particle_alternative_coordinates.velocity.y
                            / factor
                            + star_alternative_coordinates.velocity.y;
                        particle.inertial_velocity.z = particle_alternative_coordinates.velocity.z
                            / factor
                            + star_alternative_coordinates.velocity.z;
                    } else {
                        particle.inertial_velocity.x = particle_alternative_coordinates.velocity.x
                            + star_alternative_coordinates.velocity.x;
                        particle.inertial_velocity.y = particle_alternative_coordinates.velocity.y
                            + star_alternative_coordinates.velocity.y;
                        particle.inertial_velocity.z = particle_alternative_coordinates.velocity.z
                            + star_alternative_coordinates.velocity.z;
                    }
                    let factor = if self.alternative_coordinates_type == CoordinatesType::WHDS {
                        particle.mass / (m0 + particle.mass) // mf complete
                    } else {
                        particle.mass / m0
                    };
                    new_star_velocity.x -= particle_alternative_coordinates.velocity.x * factor;
                    new_star_velocity.y -= particle_alternative_coordinates.velocity.y * factor;
                    new_star_velocity.z -= particle_alternative_coordinates.velocity.z * factor;
                }
                star.inertial_velocity.x = new_star_velocity.x;
                star.inertial_velocity.y = new_star_velocity.y;
                star.inertial_velocity.z = new_star_velocity.z;
            }
        }
    }

    fn whds_and_democratic_heliocentric_to_inertial_pos(&mut self) {
        let (particles_alternative_coordinates_left, particles_alternative_coordinates_right) =
            self.particles_alternative_coordinates[..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
        if let Some((star_alternative_coordinates, particles_alternative_coordinates_right)) =
            particles_alternative_coordinates_right.split_first_mut()
        {
            let (particles_left, particles_right) = self.universe.particles
                [..self.universe.n_particles]
                .split_at_mut(self.universe.hosts.index.most_massive);
            if let Some((star, particles_right)) = particles_right.split_first_mut() {
                let mtot = star_alternative_coordinates.mass;
                let mut new_star_position = star_alternative_coordinates.position; // Copy

                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    new_star_position.x -=
                        particle_alternative_coordinates.position.x * particle.mass / mtot;
                    new_star_position.y -=
                        particle_alternative_coordinates.position.y * particle.mass / mtot;
                    new_star_position.z -=
                        particle_alternative_coordinates.position.z * particle.mass / mtot;
                }
                star.inertial_position.x = new_star_position.x;
                star.inertial_position.y = new_star_position.y;
                star.inertial_position.z = new_star_position.z;

                for (particle_alternative_coordinates, particle) in
                    particles_alternative_coordinates_left
                        .iter_mut()
                        .chain(particles_alternative_coordinates_right.iter_mut())
                        .zip(particles_left.iter_mut().chain(particles_right.iter_mut()))
                {
                    particle.inertial_position.x =
                        particle_alternative_coordinates.position.x + star.inertial_position.x;
                    particle.inertial_position.y =
                        particle_alternative_coordinates.position.y + star.inertial_position.y;
                    particle.inertial_position.z =
                        particle_alternative_coordinates.position.z + star.inertial_position.z;
                }
            }
        }
    }
}
