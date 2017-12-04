extern crate time;
use std;
use std::io::{Write, BufWriter};
use std::fs::File;
use super::Integrator;
use super::super::constants::{PI, WHFAST_NMAX_QUART, WHFAST_NMAX_NEWT, MAX_PARTICLES, G};
use super::super::particles::Universe;
use super::super::particles::IgnoreGravityTerms;
use super::super::particles::Axes;
use super::output::{write_historic_snapshot};
use rustc_serialize::json;
use bincode::rustc_serialize::{decode_from};
use bincode::SizeLimit;
use std::io::{BufReader};
use std::io::Read;
use std::path::Path;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

/// Source: Rein & Tamayo, 2015
/// WHFast is a complete reimplementation of the Wisdom-Holman integrator[1],
/// designed to speed up the algorithm and increase its accuracy (e.g., related 
/// to finite double floating-point precision, all real numbers cannot be 
/// represented exactly and floating-point precision has consequences for the 
/// numerical stability of any algorithm and the growth of numerical round-off 
/// error). It is unbiased (i.e., the errors are random and uncorrelated), and
/// has a very slow error growth. For sufficiently small timesteps, it achieves
/// Brouwer’s law (i.e., the energy error grows as time to the power of one 
/// half).
///
/// WHFast is a symplectic integrator[2], althought its symplectic nature is 
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
///     http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1612.05329
/// - Ari Silburt, Hanno Rein & Dan Tamayo
///     HERMES: a hybrid integrator for simulating close encounters and planetesimal migration
///     https://silburt.github.io/files/HERMES.pdf
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
/// [2] Source: https://en.wikipedia.org/wiki/N-body_problem#Few_bodies
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

#[derive(Debug, Copy, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub enum CoordinatesType {
    Jacobi,
    DemocraticHeliocentric,
    WHDS,
}


#[derive(Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
pub struct WHFast {
    time_step: f64,
    half_time_step: f64,
    pub universe: Universe,
    last_spin: [Axes; MAX_PARTICLES], // For spin integration with the midpoint method
    pub current_time: f64,
    current_iteration: usize,
    pub recovery_snapshot_period: f64,
    pub historic_snapshot_period: f64,
    last_recovery_snapshot_time: f64,
    last_historic_snapshot_time: f64,
    pub n_historic_snapshots: usize,
    pub hash: u64,
    /// Internal data structures below. Nothing to be changed by the user.
    particles_alternative_coordinates: [AlternativeCoordinates; MAX_PARTICLES], // Jacobi, democractic-heliocentric or WHDS
    alternative_coordinates_type: CoordinatesType,
    timestep_warning: usize ,
}

#[derive(Copy, Debug, Clone, RustcEncodable, RustcDecodable, PartialEq)]
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
        format!("{:?}", self).hash(state); 
    }
}

impl WHFast {
    pub fn new(time_step: f64, recovery_snapshot_period: f64, historic_snapshot_period: f64, universe: Universe, alternative_coordinates_type: CoordinatesType) -> WHFast {
        let particles_alternative_coordinates = [AlternativeCoordinates{mass: 0., mass_g: 0., position: Axes{x: 0., y:0., z:0.}, velocity: Axes{x: 0., y:0., z:0.}, acceleration: Axes{x: 0., y:0., z:0.}}; MAX_PARTICLES];
        let mut universe_integrator = WHFast {
                    time_step:time_step,
                    half_time_step:0.5*time_step,
                    recovery_snapshot_period:recovery_snapshot_period,
                    historic_snapshot_period:historic_snapshot_period,
                    last_recovery_snapshot_time:-1.,
                    last_historic_snapshot_time:-1.,
                    n_historic_snapshots:0,
                    hash: 0,
                    universe:universe,
                    last_spin:[Axes{x:0., y:0., z:0. }; MAX_PARTICLES],
                    current_time:0.,
                    current_iteration:0,
                    // WHFast specifics:
                    particles_alternative_coordinates: particles_alternative_coordinates,
                    alternative_coordinates_type: alternative_coordinates_type,
                    timestep_warning: 0,
                    };
        // Initialize physical values
        let current_time = 0.;
        universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
        universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
        universe_integrator
    }
    
    pub fn restore_snapshot(universe_integrator_snapshot_path: &Path, verify_integrity: bool) -> Result<WHFast, String> {
        let mut universe_integrator: WHFast;
        if universe_integrator_snapshot_path.exists() {
            // Open the path in read-only mode, returns `io::Result<File>`
            let mut snapshot_file = match File::open(&universe_integrator_snapshot_path) {
                // The `description` method of `io::Error` returns a string that
                // describes the error
                Err(why) => return Err(format!("Couldn't open {}: {}", universe_integrator_snapshot_path.display(), why)),
                Ok(file) => file,
            };


            if universe_integrator_snapshot_path.extension().unwrap() == "json" { 
                //// Deserialize using `json::decode`
                let mut json_encoded = String::new();
                match snapshot_file.read_to_string(&mut json_encoded) {
                    Err(why) => return Err(format!("Couldn't read {}: {}", universe_integrator_snapshot_path.display(), why)),
                    Ok(_) => {}
                }
                universe_integrator = json::decode(&json_encoded).unwrap();
            } else {
                let mut reader = BufReader::new(snapshot_file);
                universe_integrator = decode_from(&mut reader, SizeLimit::Infinite).unwrap();
            }
            if universe_integrator.hash == 0 {
                if verify_integrity && universe_integrator.current_time != 0. {
                    panic!("[PANIC {} UTC] File '{}' has a zeroed hash (i.e., new simulation) but a current time different from zero ({})", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display(), universe_integrator.current_time)
                }
                println!("[INFO {} UTC] Created new simulation based on '{}'.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
                // Initialize physical values
                let current_time = 0.;
                universe_integrator.universe.calculate_norm_spin(); // Needed for evolution
                universe_integrator.universe.calculate_particles_evolving_quantities(current_time); // Make sure we start with the good initial values
            } else {
                // Verify hash for this universe at this moment of time
                let mut s = DefaultHasher::new();
                let restored_hash = universe_integrator.hash;
                universe_integrator.hash = 0;
                universe_integrator.hash(&mut s);
                let computed_hash = s.finish();
                if verify_integrity && restored_hash != computed_hash {
                    panic!("[PANIC {} UTC] File '{}' seems corrupted because computed hash '{}' does not match restored hash '{}'", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display(), computed_hash, restored_hash)
                }
                universe_integrator.hash = restored_hash;
                println!("[INFO {} UTC] Restored previous simulation from '{}'.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap(), universe_integrator_snapshot_path.display());
            }
            return Ok(universe_integrator);
        } else {
            return Err(format!("File does not exist"));
        }
    }
}

impl Integrator for WHFast {

    fn iterate(&mut self, universe_history_writer: &mut BufWriter<File>, silent_mode: bool) -> Result<bool, String> {
        // Output
        let first_snapshot_trigger = self.last_historic_snapshot_time < 0.;
        let historic_snapshot_time_trigger = self.last_historic_snapshot_time + self.historic_snapshot_period <= self.current_time;
        let recovery_snapshot_time_trigger = self.last_recovery_snapshot_time + self.recovery_snapshot_period <= self.current_time;
        if first_snapshot_trigger || historic_snapshot_time_trigger {
            self.inertial_to_heliocentric_posvelacc();
            if self.universe.consider_tides {
                self.universe.calculate_denergy_dt();
            }
            write_historic_snapshot(universe_history_writer, &self.universe, self.current_time, self.time_step);
            self.last_historic_snapshot_time = self.n_historic_snapshots as f64*self.historic_snapshot_period; // Instead of self.current_time to avoid small deviations
            self.n_historic_snapshots += 1;
            let current_time_years = self.current_time/365.25;
            if ! silent_mode {
                print!("Year: {:0.0} ({:0.1e})                                              \r", current_time_years, current_time_years);
                let _ = std::io::stdout().flush();
            }
        }

        self.iterate_position_and_velocity_with_whfasthelio();
        if self.universe.evolving_particles_exist || self.universe.consider_tides || self.universe.consider_rotational_flattening {
            self.iterate_spin_with_midpoint_method();
        }

        // ---------------------------------------------------------------------
        self.current_iteration += 1;

        // Return
        if self.current_time+self.time_step > self.universe.time_limit {
            Err("reached maximum time limit.".to_string())
        } else {
            Ok(first_snapshot_trigger || recovery_snapshot_time_trigger)
        }
    }

    fn prepare_for_recovery_snapshot(&mut self, universe_history_writer: &mut BufWriter<File>) {
        self.last_recovery_snapshot_time = self.current_time;
        universe_history_writer.flush().unwrap();
        // Compute hash for this universe at this moment of time
        self.hash = 0;
        let mut s = DefaultHasher::new();
        self.hash(&mut s);
        self.hash = s.finish();
    }

}


impl WHFast {
    // WHFast integrator
    fn iterate_position_and_velocity_with_whfasthelio(&mut self) {
        let time_step = self.time_step;
        let half_time_step = self.half_time_step;

        // ---------------------------------------------------------------------
        // A 'DKD'-like integrator will do the first 'D' part.
        self.inertial_to_alternative_posvel();
        self.kepler_steps(half_time_step);
        self.jump_step(half_time_step);
        self.alternative_to_inertial_posvel();
        self.current_time += self.half_time_step;
        // ---------------------------------------------------------------------

        // Calculate accelerations.
        let ignore_gravity_terms = match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => IgnoreGravityTerms::WHFastOne,
            CoordinatesType::DemocraticHeliocentric => IgnoreGravityTerms::WHFastTwo,
            CoordinatesType::WHDS => IgnoreGravityTerms::WHFastTwo,
        };
        let ignored_gravity_terms = ignore_gravity_terms;
        self.universe.gravity_calculate_acceleration(ignore_gravity_terms);
        self.inertial_to_heliocentric_posvelacc();
        
        // Calculate spin variation and non-gravity accelerations.
        let evolution = true;
        let dspin_dt = false;
        let accelerations = true;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations, ignored_gravity_terms);
        self.heliocentric_to_inertial_posvelacc();

        // ---------------------------------------------------------------------
        // A 'DKD'-like integrator will do the 'KD' part.
        //self.inertial_to_alternative_posvel(); // not needed because position, velocities and masses have not changed
        self.interaction_step(time_step);
        self.jump_step(half_time_step);
        self.kepler_steps(half_time_step);
        self.alternative_to_inertial_posvel();
        self.current_time += self.half_time_step;
    }

    fn iterate_spin_with_midpoint_method(&mut self) {
        let time_step = self.time_step;
        let half_time_step = self.half_time_step;
        
        // ---------------------------------------------------------------------
        self.inertial_to_heliocentric_posvel();

        let ignored_gravity_terms = match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => IgnoreGravityTerms::WHFastOne,
            CoordinatesType::DemocraticHeliocentric => IgnoreGravityTerms::WHFastTwo,
            CoordinatesType::WHDS => IgnoreGravityTerms::WHFastTwo,
        };

        // Midpoint method: https://en.wikipedia.org/wiki/Midpoint_method
        let evolution = false; // Don't evolve particles at this point or it messes up the conservation of angular momentum
        let dspin_dt = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations, ignored_gravity_terms);

        self.save_current_spin();
        self.spin_step(half_time_step);
        let evolution = false; // Don't evolve particles at this point or it messes up the conservation of angular momentum
        let dspin_dt = true;
        let accelerations = false;
        self.universe.calculate_additional_effects(self.current_time, evolution, dspin_dt, accelerations, ignored_gravity_terms);
        self.restore_last_spin();
        self.spin_step(time_step);
    }

    fn save_current_spin(&mut self) {
        for (particle, spin) in self.universe.particles[..self.universe.n_particles].iter().zip(self.last_spin[..self.universe.n_particles].iter_mut()) {
            spin.x = particle.spin.x;
            spin.y = particle.spin.y;
            spin.z = particle.spin.z;
        }
    }

    fn restore_last_spin(&mut self) {
        for (particle, spin) in self.universe.particles[..self.universe.n_particles].iter_mut().zip(self.last_spin[..self.universe.n_particles].iter()) {
            particle.spin.x = spin.x;
            particle.spin.y = spin.y;
            particle.spin.z = spin.z;
        }
    }
 
    fn spin_step(&mut self, _dt: f64) {
        for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            if particle.moment_of_inertia_ratio != 1. {
                particle.spin.x = particle.moment_of_inertia_ratio * particle.spin.x + _dt * particle.dspin_dt.x;
                particle.spin.y = particle.moment_of_inertia_ratio * particle.spin.y + _dt * particle.dspin_dt.y;
                particle.spin.z = particle.moment_of_inertia_ratio * particle.spin.z + _dt * particle.dspin_dt.z;
            } else {
                particle.spin.x = particle.spin.x + _dt * particle.dspin_dt.x;
                particle.spin.y = particle.spin.y + _dt * particle.dspin_dt.y;
                particle.spin.z = particle.spin.z + _dt * particle.dspin_dt.z;
            }
            if particle.wind_factor != 0. {
                // TODO: Verify wind factor
                particle.spin.x += _dt * particle.wind_factor * particle.spin.x;
                particle.spin.y += _dt * particle.wind_factor * particle.spin.y;
                particle.spin.z += _dt * particle.wind_factor * particle.spin.z;
            }
        }
    }

    /***************************** 
     * Operators                 */
    fn jump_step(&mut self, _dt: f64){
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_jump_step(_dt),
            CoordinatesType::DemocraticHeliocentric => self.democratic_heliocentric_jump_step(_dt),
            CoordinatesType::WHDS => self.whds_jump_step(_dt),
        };
    }

    fn jacobi_jump_step(&mut self, _dt: f64){
        // Nothing to be done
    }

    fn democratic_heliocentric_jump_step(&mut self, _dt: f64){
        let m0 = self.universe.particles[0].mass;
        let mut px = 0.;
        let mut py = 0.;
        let mut pz = 0.;
        for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
            px += particle.mass* particle_alternative_coordinates.velocity.x;
            py += particle.mass* particle_alternative_coordinates.velocity.y;
            pz += particle.mass* particle_alternative_coordinates.velocity.z;
        }
        for particle_alternative_coordinates in self.particles_alternative_coordinates[1..self.universe.n_particles].iter_mut() {
            particle_alternative_coordinates.position.x += _dt * px/m0;
            particle_alternative_coordinates.position.y += _dt * py/m0;
            particle_alternative_coordinates.position.z += _dt * pz/m0;
        }
    }

    fn whds_jump_step(&mut self, _dt: f64){
        let m0 = self.universe.particles[0].mass;
        let mut px = 0.;
        let mut py = 0.;
        let mut pz = 0.;
        for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
            let factor = m0 + particle.mass; // mf numerator
            px += particle.mass* particle_alternative_coordinates.velocity.x/factor;
            py += particle.mass* particle_alternative_coordinates.velocity.y/factor;
            pz += particle.mass* particle_alternative_coordinates.velocity.z/factor;
        }
        for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter_mut().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
            let factor = m0 + particle.mass; // mf numerator
            particle_alternative_coordinates.position.x += _dt * (px - (particle.mass * particle_alternative_coordinates.velocity.x / factor));
            particle_alternative_coordinates.position.y += _dt * (py - (particle.mass * particle_alternative_coordinates.velocity.y / factor));
            particle_alternative_coordinates.position.z += _dt * (pz - (particle.mass * particle_alternative_coordinates.velocity.z / factor));
        }
    }

    fn interaction_step(&mut self, _dt: f64){
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_interaction_step(_dt),
            CoordinatesType::DemocraticHeliocentric => self.democratic_heliocentric_interaction_step(_dt),
            CoordinatesType::WHDS => self.whds_interaction_step(_dt),
        };
    }

    fn jacobi_interaction_step(&mut self, _dt: f64){
        self.inertial_to_jacobi_acc();
        let softening = 1e-12;
        let m0 = self.universe.particles[0].mass;
        if let Some((_, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
            let mut eta = m0;
            for (i, particle_alternative_coordinates) in particles_alternative_coordinates.iter_mut().enumerate() {
                eta += particle_alternative_coordinates.mass;
                particle_alternative_coordinates.velocity.x += _dt * particle_alternative_coordinates.acceleration.x;
                particle_alternative_coordinates.velocity.y += _dt * particle_alternative_coordinates.acceleration.y;
                particle_alternative_coordinates.velocity.z += _dt * particle_alternative_coordinates.acceleration.z;
                if i > 0 { // The star is already out of the loop, ignore first planet too (it should be the closer to the star)
                    let rj2i = 1./(particle_alternative_coordinates.position.x.powi(2) + particle_alternative_coordinates.position.y.powi(2) + particle_alternative_coordinates.position.z.powi(2) + softening);
                    let rji = rj2i.sqrt();
                    let rj3im = rji*rj2i*G*eta;
                    let prefac = _dt*rj3im;
                    particle_alternative_coordinates.velocity.x += prefac * particle_alternative_coordinates.position.x;
                    particle_alternative_coordinates.velocity.y += prefac * particle_alternative_coordinates.position.y;
                    particle_alternative_coordinates.velocity.z += prefac * particle_alternative_coordinates.position.z;
                }
            }
        }
    }

    fn democratic_heliocentric_interaction_step(&mut self, _dt: f64){
        // Democratic heliocentric ignores the star acceleration, it gravity is
        // always zero. But when additional forces are computed and transformed
        // to inertial coordinates, the star adquires a tiny acceleration. Thus
        // we take it into account by subtracting it to the rest of the particles
        // acceleration
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter_mut().zip(particles.iter()) {
                particle_alternative_coordinates.velocity.x += _dt*(particle.inertial_acceleration.x - star.inertial_acceleration.x);
                particle_alternative_coordinates.velocity.y += _dt*(particle.inertial_acceleration.y - star.inertial_acceleration.y);
                particle_alternative_coordinates.velocity.z += _dt*(particle.inertial_acceleration.z - star.inertial_acceleration.z);
            }
        }
    }

    fn whds_interaction_step(&mut self, _dt: f64){
        // WHDS ignores the star acceleration, it gravity is always zero. 
        // But when additional forces are computed and transformed
        // to inertial coordinates, the star adquires a tiny acceleration. Thus
        // we take it into account by subtracting it to the rest of the particles
        // acceleration
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            let m0 = star.mass;
            for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter_mut().zip(particles.iter()) {
                let factor = m0+particle.mass; // mf numerator
                particle_alternative_coordinates.velocity.x += _dt*factor*(particle.inertial_acceleration.x - star.inertial_acceleration.x)/m0;
                particle_alternative_coordinates.velocity.y += _dt*factor*(particle.inertial_acceleration.y - star.inertial_acceleration.x)/m0;
                particle_alternative_coordinates.velocity.z += _dt*factor*(particle.inertial_acceleration.z - star.inertial_acceleration.x)/m0;
            }
        }
    }


    fn kepler_steps(&mut self, time_step: f64){
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.jacobi_kepler_steps(time_step),
            CoordinatesType::DemocraticHeliocentric => self.democratic_heliocentric_kepler_steps(time_step),
            CoordinatesType::WHDS => self.whds_kepler_steps(time_step),
        };
        let star_alternative_coordinates = &mut self.particles_alternative_coordinates[0];
        star_alternative_coordinates.position.x += time_step*star_alternative_coordinates.velocity.x;
        star_alternative_coordinates.position.y += time_step*star_alternative_coordinates.velocity.y;
        star_alternative_coordinates.position.z += time_step*star_alternative_coordinates.velocity.z;
    }

    fn jacobi_kepler_steps(&mut self, time_step: f64){
        let mut star_planets_mass_g = self.universe.particles[0].mass_g;
        for i in 1..self.universe.n_particles {
            star_planets_mass_g += self.particles_alternative_coordinates[i].mass_g;
            self.kepler_individual_step(i, star_planets_mass_g, time_step);
        }
    }

    fn democratic_heliocentric_kepler_steps(&mut self, time_step: f64){
        let star_mass_g = self.universe.particles[0].mass_g;
        for i in 1..self.universe.n_particles {
            self.kepler_individual_step(i, star_mass_g, time_step);
        }
    }

    fn whds_kepler_steps(&mut self, time_step: f64){
        let star_mass_g = self.universe.particles[0].mass_g;
        for i in 1..self.universe.n_particles {
            let star_planet_mass_g = star_mass_g + self.particles_alternative_coordinates[i].mass_g;
            self.kepler_individual_step(i, star_planet_mass_g, time_step);
        }
    }

    //***************************** 
    // Keplerian motion           
    fn kepler_individual_step(&mut self, i: usize, mass_g: f64, _dt: f64){
        // Save a copy of the original position and velocities
        let p1_position;
        let p1_velocity;
        {
            let p1 = &self.particles_alternative_coordinates[i];
            p1_position = p1.position.clone();
            p1_velocity = p1.velocity.clone();
        }

        let r0 = (p1_position.x.powi(2) + p1_position.y.powi(2) + p1_position.z.powi(2)).sqrt();
        let r0i = 1./r0;
        let v2 = p1_velocity.x.powi(2) + p1_velocity.y.powi(2) + p1_velocity.z.powi(2);
        let beta = 2.*mass_g*r0i - v2;
        let eta0 = p1_position.x*p1_velocity.x + p1_position.y*p1_velocity.y + p1_position.z*p1_velocity.z;
        let zeta0 = mass_g - beta*r0;
        let mut x;
        let mut gs;
        let mut invperiod = 0.; // only used for beta>0.
        let x_per_period;
        
        if beta > 0. {
            //// Elliptic orbit
            let sqrt_beta = beta.sqrt();
            invperiod = sqrt_beta*beta / (2.*PI*mass_g);
            x_per_period = 2.*PI / sqrt_beta;
            if _dt.abs()*invperiod > 1. && self.timestep_warning == 0 {
                // Ignoring const qualifiers. This warning should not have any effect on
                // other parts of the code, nor is it vital to show it.
                println!("[WARNING {} UTC] WHFast convergence issue. Timestep is larger than at least one orbital period.", time::now_utc().strftime("%Y.%m.%d %H:%M:%S").unwrap());
                self.timestep_warning += 1;
            }
            //x = _dt*invperiod*x_per_period; // first order guess 
            let dtr0i = _dt*r0i;
            //x = dtr0i; // first order guess
            x = dtr0i * (1. - dtr0i*eta0*0.5*r0i); // second order guess
            //x = dtr0i *(1.- 0.5*dtr0i*r0i*(eta0-dtr0i*(eta0*eta0*r0i-1./3.*zeta0))); // third order guess
            //x = _dt*beta/mass_g + eta0/mass_g*(0.85*sqrt(1.+zeta0*zeta0/beta/eta0/eta0) - 1.);  // Dan's version 
        } else {
            //// Hyperbolic orbit
            x = 0.; // Initial guess 
            x_per_period = std::f64::NAN; // only used for beta>0. nan triggers Newton's method for beta<0
        }


        let mut converged = 0;
        let mut old_x = x;

        //// Do one Newton step
        gs = WHFast::stiefel_gs3(beta, x);
        let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
        let mut ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        x  = ri*(x*eta0_gs1_zeta0_gs2-eta0*gs[2]-zeta0*gs[3]+_dt);

        // Choose solver depending on estimated step size
        // Note, for hyperbolic orbits this uses Newton's method.
        if !x_per_period.is_nan() && (x-old_x).abs() > 0.01*x_per_period {
            // Quartic solver
            // Linear initial guess
            x = beta*_dt/mass_g;
            let mut prev_x = [0.; WHFAST_NMAX_QUART+1];
            'outer: for n_lag in 1..WHFAST_NMAX_QUART {
                gs = WHFast::stiefel_gs3(beta, x);
                let f = r0*x + eta0*gs[2] + zeta0*gs[3] - _dt;
                let fp = r0 + eta0*gs[1] + zeta0*gs[2];
                let fpp = eta0*gs[0] + zeta0*gs[1];
                let denom = fp + (16.*fp*fp - 20.*f*fpp).abs().sqrt();
                x = (x*denom - 5.*f)/denom;
                'inner: for i in 1..n_lag {
                    if x == prev_x[i] {
                        // Converged. Exit.
                        //n_lag = WHFAST_NMAX_QUART;
                        converged = 1;
                        break 'outer;
                    }
                }
                prev_x[n_lag] = x;
            }
            let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
            ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        } else {
            // Newton's method
            let mut old_x2;
            for _ in 1..WHFAST_NMAX_NEWT {
                old_x2 = old_x;
                old_x = x;
                gs = WHFast::stiefel_gs3(beta, x);
                let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
                ri = 1./(r0 + eta0_gs1_zeta0_gs2);
                x  = ri*(x*eta0_gs1_zeta0_gs2-eta0*gs[2]-zeta0*gs[3]+_dt);
                
                if x==old_x || x==old_x2 {
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
                x_min = x_per_period * (_dt*invperiod).floor();
                x_max = x_min + x_per_period;
            } else {
                //Hyperbolic
                let h2 = r0*r0*v2-eta0*eta0;
                let q = h2/mass_g/(1.+(1.-h2*beta/(mass_g*mass_g)).sqrt());
                let vq = (h2).sqrt()/q;
                x_min = 1./(vq+r0/_dt);
                x_max = _dt/q;
            }
            x = (x_max + x_min)/2.;
            loop {
                gs = WHFast::stiefel_gs3(beta, x);
                let s   = r0*x + eta0*gs[2] + zeta0*gs[3]-_dt;
                if s >= 0. {
                    x_max = x;
                } else {
                    x_min = x;
                }
                x = (x_max + x_min)/2.;

                if (x_max-x_min).abs()/x_max <= 1e-15 {
                    break;
                }
            }
            let eta0_gs1_zeta0_gs2 = eta0*gs[1] + zeta0*gs[2];
            ri = 1./(r0 + eta0_gs1_zeta0_gs2);
        }
        if ri.is_nan() {
            // Exception for (almost) straight line motion in hyperbolic case
            ri = 0.;
            gs[1] = 0.;
            gs[2] = 0.;
            gs[3] = 0.;
        }

        // Note: These are not the traditional f and g functions.
        let f = -mass_g*gs[2]*r0i;
        let g = _dt - mass_g*gs[3];
        let fd = -mass_g*gs[1]*r0i*ri; 
        let gd = -mass_g*gs[2]*ri; 
            
        let p_j = &mut self.particles_alternative_coordinates[i];
        p_j.position.x += f*p1_position.x + g*p1_velocity.x;
        p_j.position.y += f*p1_position.y + g*p1_velocity.y;
        p_j.position.z += f*p1_position.z + g*p1_velocity.z;
        
        // WARNING: p1_position used below should not be modified by the previous block (keep p_j
        // independent)
        p_j.velocity.x += fd*p1_position.x + gd*p1_velocity.x;
        p_j.velocity.y += fd*p1_position.y + gd*p1_velocity.y;
        p_j.velocity.z += fd*p1_position.z + gd*p1_velocity.z;
    }

    fn stiefel_gs3(beta: f64, x: f64) -> [f64; 6] {
        let x2 = x.powi(2);
        let mut gs = WHFast::stumpff_cs3(beta*x2);
        gs[1] *= x; 
        gs[2] *= x2; 
        gs[3] *= x2*x;
        gs
    }

    fn stumpff_cs3(z: f64) -> [f64; 6] {
        // Fast inverse factorial lookup table
        let invfactorial: [f64; 35] = [1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040., 1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800., 1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000., 1./6402373705728000., 1./121645100408832000., 1./2432902008176640000., 1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000., 1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000., 1./10888869450418352160768000000., 1./304888344611713860501504000000., 1./8841761993739701954543616000000., 1./265252859812191058636308480000000., 1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000., 1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.];
        let mut z = z;
        let mut n = 0;
        while z.abs() > 0.1 {
            z = z/4.;
            n += 1;
        }
        let mut cs = [0.; 6];
        let nmax = 13;
        let mut c_odd = invfactorial[nmax];
        let mut c_even = invfactorial[nmax-1];

        let mut np = nmax-2;
        while np >= 3 {
            c_odd  = invfactorial[np]    - z *c_odd;
            c_even = invfactorial[np-1]  - z *c_even;
            np -= 2;
        }
        cs[3] = c_odd;
        cs[2] = c_even;
        cs[1] = invfactorial[1]  - z *c_odd;
        cs[0] = invfactorial[0]  - z *c_even;
        while n > 0 {
            cs[3] = (cs[2]+cs[0]*cs[3])*0.25;
            cs[2] = cs[1]*cs[1]*0.5;
            cs[1] = cs[0]*cs[1];
            cs[0] = 2.*cs[0]*cs[0]-1.;
            n -= 1;
        }
        cs
    }

    //***************************** 
    // Coordinate transformations 
    //***************************** 
    fn inertial_to_alternative_posvel(&mut self){
        match self.alternative_coordinates_type {
            CoordinatesType::Jacobi => self.inertial_to_jacobi_posvel(),
            CoordinatesType::DemocraticHeliocentric => self.inertial_to_democratic_heliocentric_posvel(),
            CoordinatesType::WHDS => self.inertial_to_whds_posvel(),
        };
    }

    fn inertial_to_jacobi_posvel(&mut self){
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            let m0 = star.mass;
            let mut eta = m0;
            let mut s_x = eta * star.inertial_position.x;
            let mut s_y = eta * star.inertial_position.y;
            let mut s_z = eta * star.inertial_position.z;
            let mut s_vx = eta * star.inertial_velocity.x;
            let mut s_vy = eta * star.inertial_velocity.y;
            let mut s_vz = eta * star.inertial_velocity.z;
            if let Some((star_alternative_coordinates, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
                for (particle_alternative_coordinates, particle) in particles_alternative_coordinates.iter_mut().zip(particles.iter()) {
                    let ei = 1./eta;
                    eta += particle.mass;
                    let pme = eta*ei;
                    particle_alternative_coordinates.mass = particle.mass;
                    particle_alternative_coordinates.mass_g = particle.mass_g;
                    particle_alternative_coordinates.position.x = particle.inertial_position.x - s_x*ei;
                    particle_alternative_coordinates.position.y = particle.inertial_position.y - s_y*ei;
                    particle_alternative_coordinates.position.z = particle.inertial_position.z - s_z*ei;
                    particle_alternative_coordinates.velocity.x = particle.inertial_velocity.x - s_vx*ei;
                    particle_alternative_coordinates.velocity.y = particle.inertial_velocity.y - s_vy*ei;
                    particle_alternative_coordinates.velocity.z = particle.inertial_velocity.z - s_vz*ei;
                    s_x  = s_x  * pme + particle.mass*particle_alternative_coordinates.position.x ;
                    s_y  = s_y  * pme + particle.mass*particle_alternative_coordinates.position.y ;
                    s_z  = s_z  * pme + particle.mass*particle_alternative_coordinates.position.z ;
                    s_vx = s_vx * pme + particle.mass*particle_alternative_coordinates.velocity.x;
                    s_vy = s_vy * pme + particle.mass*particle_alternative_coordinates.velocity.y;
                    s_vz = s_vz * pme + particle.mass*particle_alternative_coordinates.velocity.z;
                }
                let mtot = eta;
                let mtot_i = 1./mtot;
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

    fn inertial_to_jacobi_acc(&mut self){
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            if let Some((star_alternative_coordinates, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
                let mut eta = star.mass;
                let mut s_ax = eta * star.inertial_acceleration.x;
                let mut s_ay = eta * star.inertial_acceleration.y;
                let mut s_az = eta * star.inertial_acceleration.z;
                for (particle_alternative_coordinates, particle) in particles_alternative_coordinates.iter_mut().zip(particles.iter()) {
                    let ei = 1./eta;
                    eta += particle.mass;
                    let pme = eta*ei;
                    particle_alternative_coordinates.acceleration.x = particle.inertial_acceleration.x - s_ax*ei;
                    particle_alternative_coordinates.acceleration.y = particle.inertial_acceleration.y - s_ay*ei;
                    particle_alternative_coordinates.acceleration.z = particle.inertial_acceleration.z - s_az*ei;
                    s_ax = s_ax * pme + particle.mass*particle_alternative_coordinates.acceleration.x;
                    s_ay = s_ay * pme + particle.mass*particle_alternative_coordinates.acceleration.y;
                    s_az = s_az * pme + particle.mass*particle_alternative_coordinates.acceleration.z;
                }
                let mtot = eta;
                let mtot_i = 1./mtot;
                star_alternative_coordinates.acceleration.x = s_ax * mtot_i;
                star_alternative_coordinates.acceleration.y = s_ay * mtot_i;
                star_alternative_coordinates.acceleration.z = s_az * mtot_i;
            }
        }
    }

    fn inertial_to_democratic_heliocentric_posvel(&mut self){
        self.inertial_to_whds_and_democratic_heliocentric_posvel();
    }

    fn inertial_to_whds_posvel(&mut self){
        self.inertial_to_whds_and_democratic_heliocentric_posvel();
    }

    fn inertial_to_whds_and_democratic_heliocentric_posvel(&mut self){
        {
            let star_alternative_coordinates = &mut self.particles_alternative_coordinates[0];
            star_alternative_coordinates.position.x = 0.;
            star_alternative_coordinates.position.y = 0.;
            star_alternative_coordinates.position.z = 0.;
            star_alternative_coordinates.velocity.x = 0.;
            star_alternative_coordinates.velocity.y = 0.;
            star_alternative_coordinates.velocity.z = 0.;
            star_alternative_coordinates.mass = 0.;
            star_alternative_coordinates.mass_g = 0.;
            for particle in self.universe.particles[0..self.universe.n_particles].iter() {
                star_alternative_coordinates.position.x += particle.inertial_position.x*particle.mass;
                star_alternative_coordinates.position.y += particle.inertial_position.y*particle.mass;
                star_alternative_coordinates.position.z += particle.inertial_position.z*particle.mass;
                star_alternative_coordinates.velocity.x += particle.inertial_velocity.x*particle.mass;
                star_alternative_coordinates.velocity.y += particle.inertial_velocity.y*particle.mass;
                star_alternative_coordinates.velocity.z += particle.inertial_velocity.z*particle.mass;
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
        }

        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            if let Some((star_alternative_coordinates, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
                for (particle_alternative_coordinates, particle) in particles_alternative_coordinates.iter_mut().zip(particles.iter()) {
                    particle_alternative_coordinates.position.x = particle.inertial_position.x - star.inertial_position.x ; 
                    particle_alternative_coordinates.position.y = particle.inertial_position.y - star.inertial_position.y ;
                    particle_alternative_coordinates.position.z = particle.inertial_position.z - star.inertial_position.z ;
                    particle_alternative_coordinates.velocity.x = particle.inertial_velocity.x - star_alternative_coordinates.velocity.x;
                    particle_alternative_coordinates.velocity.y = particle.inertial_velocity.y - star_alternative_coordinates.velocity.y;
                    particle_alternative_coordinates.velocity.z = particle.inertial_velocity.z - star_alternative_coordinates.velocity.z;
                    particle_alternative_coordinates.mass = particle.mass;
                    particle_alternative_coordinates.mass_g = particle.mass_g;
                    if self.alternative_coordinates_type == CoordinatesType::WHDS {
                        let factor = (star.mass+particle.mass)/star.mass; // mf complete
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
            CoordinatesType::DemocraticHeliocentric => self.democratic_heliocentric_to_inertial_posvel(),
            CoordinatesType::WHDS => self.whds_to_inertial_posvel(),
        };
    }

    fn jacobi_to_inertial_posvel(&mut self) {
        if let Some((star_alternative_coordinates, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
            let m0 = star_alternative_coordinates.mass;
            let mut eta = m0;
            let mut s_x = eta * star_alternative_coordinates.position.x;
            let mut s_y = eta * star_alternative_coordinates.position.y;
            let mut s_z = eta * star_alternative_coordinates.position.z;
            let mut s_vx = eta * star_alternative_coordinates.velocity.x;
            let mut s_vy = eta * star_alternative_coordinates.velocity.y;
            let mut s_vz = eta * star_alternative_coordinates.velocity.z;
            if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
                for (particle, particle_alternative_coordinates) in particles.iter_mut().rev().zip(particles_alternative_coordinates.iter().rev()) {
                    let ei = 1./eta;
                    s_x  = (s_x  - particle.mass*particle_alternative_coordinates.position.x) * ei;
                    s_y  = (s_y  - particle.mass*particle_alternative_coordinates.position.y) * ei;
                    s_z  = (s_z  - particle.mass*particle_alternative_coordinates.position.z) * ei;
                    s_vx = (s_vx - particle.mass*particle_alternative_coordinates.velocity.x) * ei;
                    s_vy = (s_vy - particle.mass*particle_alternative_coordinates.velocity.y) * ei;
                    s_vz = (s_vz - particle.mass*particle_alternative_coordinates.velocity.z) * ei;
                    particle.inertial_position.x = particle_alternative_coordinates.position.x + s_x;
                    particle.inertial_position.y = particle_alternative_coordinates.position.y + s_y;
                    particle.inertial_position.z = particle_alternative_coordinates.position.z + s_z;
                    particle.inertial_velocity.x = particle_alternative_coordinates.velocity.x + s_vx;
                    particle.inertial_velocity.y = particle_alternative_coordinates.velocity.y + s_vy;
                    particle.inertial_velocity.z = particle_alternative_coordinates.velocity.z + s_vz;
                    eta -= particle.mass;
                    s_x  *= eta;
                    s_y  *= eta;
                    s_z  *= eta;
                    s_vx *= eta;
                    s_vy *= eta;
                    s_vz *= eta;
                }
                let mtot = eta;
                let mtot_i = 1./mtot;
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
        let m0 = self.universe.particles[0].mass;
        if let Some((star_alternative_coordinates, particles_alternative_coordinates)) = self.particles_alternative_coordinates[..self.universe.n_particles].split_first_mut() {
            for (particle_alternative_coordinates, particle) in particles_alternative_coordinates.iter().zip(self.universe.particles[1..self.universe.n_particles].iter_mut()) {
                if self.alternative_coordinates_type == CoordinatesType::WHDS {
                    let factor = (m0+particle.mass)/m0; // mf complete
                    particle.inertial_velocity.x = particle_alternative_coordinates.velocity.x/factor + star_alternative_coordinates.velocity.x;
                    particle.inertial_velocity.y = particle_alternative_coordinates.velocity.y/factor + star_alternative_coordinates.velocity.y;
                    particle.inertial_velocity.z = particle_alternative_coordinates.velocity.z/factor + star_alternative_coordinates.velocity.z;
                } else {
                    particle.inertial_velocity.x = particle_alternative_coordinates.velocity.x + star_alternative_coordinates.velocity.x;
                    particle.inertial_velocity.y = particle_alternative_coordinates.velocity.y + star_alternative_coordinates.velocity.y;
                    particle.inertial_velocity.z = particle_alternative_coordinates.velocity.z + star_alternative_coordinates.velocity.z;
                }
            }
        }

        let mut new_star_velocity = self.particles_alternative_coordinates[0].velocity.clone();
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for (particle, particle_alternative_coordinates) in particles.iter().zip(self.particles_alternative_coordinates[1..self.universe.n_particles].iter()) {
                let factor;
                if self.alternative_coordinates_type == CoordinatesType::WHDS {
                    factor = particle.mass/(m0+particle.mass) // mf complete
                } else {
                    factor = particle.mass/m0;
                }
                new_star_velocity.x -= particle_alternative_coordinates.velocity.x*factor;
                new_star_velocity.y -= particle_alternative_coordinates.velocity.y*factor;
                new_star_velocity.z -= particle_alternative_coordinates.velocity.z*factor;
            }
            star.inertial_velocity.x = new_star_velocity.x;
            star.inertial_velocity.y = new_star_velocity.y;
            star.inertial_velocity.z = new_star_velocity.z;
        }
    }

    fn whds_and_democratic_heliocentric_to_inertial_pos(&mut self) {
        let mtot = self.particles_alternative_coordinates[0].mass;
        {
            let mut new_star_position = self.particles_alternative_coordinates[0].position; // Copy
            for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter().zip(self.universe.particles[1..self.universe.n_particles].iter()) {
                new_star_position.x  -= particle_alternative_coordinates.position.x*particle.mass/mtot;
                new_star_position.y  -= particle_alternative_coordinates.position.y*particle.mass/mtot;
                new_star_position.z  -= particle_alternative_coordinates.position.z*particle.mass/mtot;
            }
            let star = &mut self.universe.particles[0];
            star.inertial_position.x  = new_star_position.x;
            star.inertial_position.y  = new_star_position.y;
            star.inertial_position.z  = new_star_position.z;
        }
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for (particle_alternative_coordinates, particle) in self.particles_alternative_coordinates[1..self.universe.n_particles].iter().zip(particles.iter_mut()) {
                particle.inertial_position.x = particle_alternative_coordinates.position.x + star.inertial_position.x;
                particle.inertial_position.y = particle_alternative_coordinates.position.y + star.inertial_position.y;
                particle.inertial_position.z = particle_alternative_coordinates.position.z + star.inertial_position.z;
            }
        }

    }

    fn inertial_to_heliocentric_posvelacc(&mut self) {
        self.inertial_to_heliocentric_posvel();
        self.inertial_to_heliocentric_acc();
        //// Equivalent formulation to substractinga stellar inertial position/velocity/acceleration:
        ////
        //// heliocentric_position = inertial_position + (1/star.mass) * sum(particle.mass * particle.position)
        //// heliocentric_velocity = inertial_velocity + (1/star.mass) * sum(particle.mass * particle.velocity)
        //// heliocentric_acceleration = inertial_acceleration + (1/star.mass) * sum(particle.mass * particle.acceleration)
        ////
        //let mut reference_position = Axes{x:0., y:0., z:0.};
        //let mut reference_velocity = Axes{x:0., y:0., z:0.};
        //let mut reference_acceleration = Axes{x:0., y:0., z:0.};
        //let star_mass = self.universe.particles[0].mass;

        //for particle in self.universe.particles[..self.universe.n_particles].iter() {
            //reference_position.x      += particle.position.x*particle.mass;
            //reference_position.y      += particle.position.y*particle.mass;
            //reference_position.z      += particle.position.z*particle.mass;
            //reference_velocity.x      += particle.velocity.x*particle.mass;
            //reference_velocity.y      += particle.velocity.y*particle.mass;
            //reference_velocity.z      += particle.velocity.z*particle.mass;
            //reference_acceleration.x  += particle.acceleration.x*particle.mass;
            //reference_acceleration.y  += particle.acceleration.y*particle.mass;
            //reference_acceleration.z  += particle.acceleration.z*particle.mass;
        //}
        //reference_position.x /= star_mass;
        //reference_position.y /= star_mass;
        //reference_position.z /= star_mass;
        //reference_velocity.x /= star_mass;
        //reference_velocity.y /= star_mass;
        //reference_velocity.z /= star_mass;
        //reference_acceleration.x /= star_mass;
        //reference_acceleration.y /= star_mass;
        //reference_acceleration.z /= star_mass;

        //for particle in self.universe.particles[..self.universe.n_particles].iter_mut() {
            //println!("++ before position {:?}", particle.position);
            //particle.position.x = particle.inertial_position.x + reference_position.x;
            //particle.position.y = particle.inertial_position.y + reference_position.y;
            //particle.position.z = particle.inertial_position.z + reference_position.z;
            //println!("++ after position {:?}", particle.position);
            //println!("++ before velocity {:?}", particle.velocity);
            //particle.velocity.x = particle.inertial_velocity.x + reference_velocity.x;
            //particle.velocity.y = particle.inertial_velocity.y + reference_velocity.y;
            //particle.velocity.z = particle.inertial_velocity.z + reference_velocity.z;
            //println!("++ after velocity {:?}", particle.velocity);
            //println!("++ before acceleration {:?}", particle.velocity);
            //particle.acceleration.x = particle.inertial_acceleration.x + reference_acceleration.x;
            //particle.acceleration.y = particle.inertial_acceleration.y + reference_acceleration.y;
            //particle.acceleration.z = particle.inertial_acceleration.z + reference_acceleration.z;
            //println!("++ after acceleration {:?}", particle.velocity);
        //}
        //println!("-----------------");
    }

    fn inertial_to_heliocentric_posvel(&mut self) {
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.position.x = particle.inertial_position.x - star.inertial_position.x;
                particle.position.y = particle.inertial_position.y - star.inertial_position.y;
                particle.position.z = particle.inertial_position.z - star.inertial_position.z;
                particle.velocity.x = particle.inertial_velocity.x - star.inertial_velocity.x;
                particle.velocity.y = particle.inertial_velocity.y - star.inertial_velocity.y;
                particle.velocity.z = particle.inertial_velocity.z - star.inertial_velocity.z;
            }
            star.position.x = 0.;
            star.position.y = 0.;
            star.position.z = 0.;
            star.velocity.x = 0.;
            star.velocity.y = 0.;
            star.velocity.z = 0.;
        }
    }

    fn inertial_to_heliocentric_acc(&mut self) {
        if let Some((star, particles)) = self.universe.particles[..self.universe.n_particles].split_first_mut() {
            for particle in particles.iter_mut() {
                particle.acceleration.x = particle.inertial_acceleration.x - star.inertial_acceleration.x;
                particle.acceleration.y = particle.inertial_acceleration.y - star.inertial_acceleration.y;
                particle.acceleration.z = particle.inertial_acceleration.z - star.inertial_acceleration.z;
            }
            star.acceleration.x = 0.;
            star.acceleration.y = 0.;
            star.acceleration.z = 0.;
        }
    }


    fn heliocentric_to_inertial_posvelacc(&mut self) {
        // Formulation:
        //
        // inertial_position = heliocentric_position - (1/total_mass) * sum(particle.mass * particle.position)
        // inertial_velocity = heliocentric_velocity - (1/total_mass) * sum(particle.mass * particle.velocity)
        // inertial_acceleration = heliocentric_acceleration - (1/total_mass) * sum(particle.mass * particle.acceleration)
        //
        // Compute center of mass
        let mut reference_position = Axes{x:0., y:0., z:0.};
        let mut reference_velocity = Axes{x:0., y:0., z:0.};
        let mut reference_acceleration = Axes{x:0., y:0., z:0.};
        let mut total_mass = 0.;

        for particle in self.universe.particles[..self.universe.n_particles].iter() {
            reference_position.x      += particle.position.x*particle.mass;
            reference_position.y      += particle.position.y*particle.mass;
            reference_position.z      += particle.position.z*particle.mass;
            reference_velocity.x      += particle.velocity.x*particle.mass;
            reference_velocity.y      += particle.velocity.y*particle.mass;
            reference_velocity.z      += particle.velocity.z*particle.mass;
            reference_acceleration.x  += particle.acceleration.x*particle.mass;
            reference_acceleration.y  += particle.acceleration.y*particle.mass;
            reference_acceleration.z  += particle.acceleration.z*particle.mass;
            total_mass += particle.mass;
        }
        reference_position.x /= total_mass;
        reference_position.y /= total_mass;
        reference_position.z /= total_mass;
        reference_velocity.x /= total_mass;
        reference_velocity.y /= total_mass;
        reference_velocity.z /= total_mass;
        reference_acceleration.x /= total_mass;
        reference_acceleration.y /= total_mass;
        reference_acceleration.z /= total_mass;

        // Convert heliocentric to inertial
        for particle in self.universe.particles[0..self.universe.n_particles].iter_mut() {
            particle.inertial_position.x = particle.position.x - reference_position.x;
            particle.inertial_position.y = particle.position.y - reference_position.y;
            particle.inertial_position.z = particle.position.z - reference_position.z;
            particle.inertial_velocity.x = particle.velocity.x - reference_velocity.x;
            particle.inertial_velocity.y = particle.velocity.y - reference_velocity.y;
            particle.inertial_velocity.z = particle.velocity.z - reference_velocity.z;
            particle.inertial_acceleration.x = particle.acceleration.x - reference_acceleration.x;
            particle.inertial_acceleration.y = particle.acceleration.y - reference_acceleration.y;
            particle.inertial_acceleration.z = particle.acceleration.z - reference_acceleration.z;
        }
    }

}
