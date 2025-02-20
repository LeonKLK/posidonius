use super::super::constants::{G, MAX_DISTANCE_2, MAX_PARTICLES};
use super::super::effects::{
    disk, evolution, general_relativity, rotational_flattening, tides, wind,
};
use super::super::{DiskEffect, RotationalFlatteningEffect, TidesEffect, WindEffect};
use super::super::{EvolutionType, Evolver};
use super::super::{GeneralRelativityEffect, GeneralRelativityImplementation};
use super::Axes;
use super::Particle;
use super::common;
use itertools::izip;
use serde::{Deserialize, Serialize};
use serde_big_array::BigArray;
use std::collections::HashMap;
use time;
use time::{OffsetDateTime, format_description};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HostIndices {
    pub most_massive: usize,
    tides: usize,                 // Particle that is the main one for tidal effects
    rotational_flattening: usize, // Particle that is the main one for rotational flattenning effects
    general_relativity: usize,    // Central particle for general relativity effects
    disk: usize,                  // Particle with disk
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct HostMostMassive {
    all: bool, // True if all are the most massive (all have the same index)
    general_relativity: bool,
    tides: bool,
    rotational_flattening: bool,
    disk: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Hosts {
    pub index: HostIndices,        // Position in the particles array
    most_massive: HostMostMassive, // Optimization: Is a particular host also the most massive?
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ConsiderEffects {
    pub tides: bool,
    pub rotational_flattening: bool,
    pub general_relativity: bool,
    pub disk: bool,
    pub wind: bool,
    pub evolution: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Universe {
    pub initial_time: f64,
    pub time_limit: f64,
    #[serde(with = "BigArray")]
    pub particles: [Particle; MAX_PARTICLES],
    pub particles_evolvers: Vec<Evolver>,
    pub n_particles: usize,
    pub consider_effects: ConsiderEffects,
    pub general_relativity_implementation: GeneralRelativityImplementation, // Optimization: fast access to GR implementation for integrators
    pub hosts: Hosts,
    pair_dependent_scaled_dissipation_factor: HashMap<usize, f64>, // Central body specific
    #[serde(with = "BigArray")]
    roche_radiuses: [f64; MAX_PARTICLES * MAX_PARTICLES],
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum IgnoreGravityTerms {
    None,
    WHFastOne,
    WHFastTwo,
}

impl Universe {
    pub fn new(
        initial_time: f64,
        time_limit: f64,
        mut particles: Vec<Particle>,
        mut consider_effects: ConsiderEffects,
    ) -> Universe {
        disable_unnecessary_effects(&mut consider_effects, &particles);
        check_effects_vs_central_and_orbiting(&particles, &consider_effects);
        let hosts = find_indices(&particles, &consider_effects);

        // Initialize general relativity factor
        let mut general_relativity_implementation = GeneralRelativityImplementation::Disabled;
        if consider_effects.general_relativity {
            let (particles_left, particles_right) =
                particles.split_at_mut(hosts.index.general_relativity);
            if let Some((general_relativity_host_particle, particles_right)) =
                particles_right.split_first_mut()
            {
                if let GeneralRelativityEffect::CentralBody(implementation) =
                    general_relativity_host_particle.general_relativity.effect
                {
                    if implementation != GeneralRelativityImplementation::Disabled {
                        general_relativity_implementation = implementation;
                        let local_copy_star_mass_g = general_relativity_host_particle.mass_g;
                        for particle in particles_left.iter_mut().chain(particles_right.iter_mut())
                        {
                            particle.general_relativity.parameters.internal.factor =
                                local_copy_star_mass_g * particle.mass_g
                                    / (local_copy_star_mass_g + particle.mass_g).powi(2);
                        }
                    }
                }
            }
        }

        // Re-compute inertial positions/velocities (if the user introduced heliocentric
        // positions/velocities, they are transformed to barycentric)
        let (center_of_mass_position, center_of_mass_velocity) =
            calculate_center_of_mass(&particles);
        for particle in &mut particles {
            particle.inertial_position.x =
                particle.heliocentric_position.x - center_of_mass_position.x;
            particle.inertial_position.y =
                particle.heliocentric_position.y - center_of_mass_position.y;
            particle.inertial_position.z =
                particle.heliocentric_position.z - center_of_mass_position.z;
            particle.inertial_velocity.x =
                particle.heliocentric_velocity.x - center_of_mass_velocity.x;
            particle.inertial_velocity.y =
                particle.heliocentric_velocity.y - center_of_mass_velocity.y;
            particle.inertial_velocity.z =
                particle.heliocentric_velocity.z - center_of_mass_velocity.z;
        }

        // OPTIMIZATION: Transform vector to array
        // - Arrays are stored in the stack which is faster than the heap (where vectors are allocated)
        // - The array should have a fixed size, thus it should always be equal or greater to the vector passed
        // - The array elements to be considered will be limited by the n_particles value
        let n_particles = particles.len();
        assert!(
            (n_particles <= MAX_PARTICLES),
            "Only {MAX_PARTICLES} bodies are allowed, you need to increase the MAX_PARTICLE constant."
        );
        let mut transformed_particles = [Particle::new_dummy(); MAX_PARTICLES];
        let mut particles_evolvers: Vec<Evolver> = Vec::with_capacity(n_particles);
        for i in 0..n_particles {
            transformed_particles[i] = particles[i];
            transformed_particles[i].id = i;
            particles_evolvers.push(Evolver::new(
                transformed_particles[i].evolution,
                initial_time,
                time_limit,
            ));
        }
        for _i in n_particles..MAX_PARTICLES {
            // For dummy particles
            particles_evolvers.push(Evolver::new(
                EvolutionType::NonEvolving,
                initial_time,
                time_limit,
            ));
        }

        if consider_effects.evolution {
            // Check if moment of inertia and angular momentum needs to be recomputed to match the
            // proper state of the selected evolving body model (i.e., initial radius and radius of
            // gyration are good given the initial time)
            for (particle, evolver) in izip!(&mut transformed_particles, &mut particles_evolvers) {
                if EvolutionType::NonEvolving == particle.evolution {
                    continue;
                }
                let mut update_angular_momentum = false;
                let current_time = 0.;
                let new_radius = evolver.radius(current_time, particle.radius);
                if abs!(new_radius - particle.radius) > 1e-6 {
                    println!(
                        "[WARNING {} UTC] Changed radius value from '{:.6}' to '{:.6}' to match expected state following the selected evolving body model",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap(),
                        particle.radius,
                        new_radius
                    );
                    particle.radius = new_radius;
                    update_angular_momentum = true;
                }
                let new_radius_of_gyration_2 =
                    evolver.radius_of_gyration_2(current_time, particle.radius_of_gyration_2);
                if abs!(new_radius_of_gyration_2 - particle.radius_of_gyration_2) > 1e-6 {
                    println!(
                        "[WARNING {} UTC] Changed radius of gyration value from '{:.6}' to '{:.6}' to match expected state following the selected evolving body model",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap(),
                        sqrt!(particle.radius_of_gyration_2),
                        sqrt!(new_radius_of_gyration_2)
                    );
                    particle.radius_of_gyration_2 = new_radius_of_gyration_2;
                    update_angular_momentum = true;
                }
                if update_angular_momentum {
                    println!(
                        "[WARNING {} UTC] Recomputed moment of inertia and angular momentum to match expected state following the selected evolving body model",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap()
                    );
                    particle.moment_of_inertia =
                        particle.mass * particle.radius_of_gyration_2 * particle.radius.powi(2);
                    particle.angular_momentum.x = particle.spin.x * particle.moment_of_inertia;
                    particle.angular_momentum.y = particle.spin.y * particle.moment_of_inertia;
                    particle.angular_momentum.z = particle.spin.z * particle.moment_of_inertia;
                }
            }
        }

        let roche_radiuses: [f64; MAX_PARTICLES * MAX_PARTICLES] =
            [0.; MAX_PARTICLES * MAX_PARTICLES];

        Universe {
            initial_time,
            time_limit,
            particles: transformed_particles,
            particles_evolvers,
            n_particles,
            consider_effects,
            general_relativity_implementation,
            hosts,
            pair_dependent_scaled_dissipation_factor: HashMap::new(),
            roche_radiuses,
        }
    }

    pub fn calculate_roche_radiuses(&mut self) {
        let (particles, _) = self.particles.split_at_mut(self.n_particles);
        let (roche_radiuses, _) = self
            .roche_radiuses
            .split_at_mut(self.n_particles * self.n_particles);
        for (i, (particle_a, roche_radiuses)) in particles
            .iter()
            .zip(roche_radiuses.chunks_mut(self.n_particles))
            .enumerate()
        {
            for (j, (particle_b, roche_radius)) in
                particles.iter().zip(roche_radiuses.iter_mut()).enumerate()
            {
                // Roche radius calculation
                if i == j {
                    continue;
                }
                // Faber et al, 2005; Pacynski, 1971
                if particle_a.mass > particle_b.mass {
                    // particle a is the most massive of both
                    *roche_radius = (particle_b.radius / 0.462)
                        * (particle_a.mass / particle_b.mass).powf(1. / 3.);
                } else {
                    // particle b is the most massive of both
                    *roche_radius = (particle_a.radius / 0.462)
                        * (particle_b.mass / particle_a.mass).powf(1. / 3.);
                }
            }
        }
    }

    pub fn gravity_calculate_acceleration(&mut self, ignore_terms: IgnoreGravityTerms) {
        let (particles, _) = self.particles.split_at_mut(self.n_particles);
        //let (roche_radiuses, _) = self.roche_radiuses.split_at_mut(self.n_particles*);
        let (roche_radiuses, _) = self
            .roche_radiuses
            .split_at_mut(self.n_particles * self.n_particles);

        let mut newtonian_inertial_accelerations = [Axes::new(); MAX_PARTICLES];
        let (newtonian_inertial_accelerations, _) =
            newtonian_inertial_accelerations.split_at_mut(self.n_particles);

        //// For compensated summation:
        //let mut newtonian_inertial_acceleration_errors = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        //let (newtonian_inertial_acceleration_errors, _) = newtonian_inertial_acceleration_errors.split_at_mut(self.n_particles);

        //// For compensated summation:
        //for (i, (particle_a, (newtonian_inertial_acceleration, (newtonian_inertial_acceleration_error, roche_radiuses)))) in particles.iter().zip(newtonian_inertial_accelerations.iter_mut().zip(newtonian_inertial_acceleration_errors.iter_mut().zip(roche_radiuses.iter()))).enumerate() {
        for (i, (particle_a, (newtonian_inertial_acceleration, roche_radiuses))) in particles
            .iter()
            .zip(
                newtonian_inertial_accelerations
                    .iter_mut()
                    .zip(roche_radiuses.chunks_mut(self.n_particles)),
            )
            .enumerate()
        {
            for (j, (particle_b, roche_radius)) in
                particles.iter().zip(roche_radiuses.iter()).enumerate()
            {
                if i == j {
                    continue;
                }

                //// Check collisions and ejections //////////////////////////////////
                let dx = particle_a.inertial_position.x - particle_b.inertial_position.x;
                let dy = particle_a.inertial_position.y - particle_b.inertial_position.y;
                let dz = particle_a.inertial_position.z - particle_b.inertial_position.z;
                let distance_2 = dx * dx + dy * dy + dz * dz;
                // Do not check twice the same pair of particles:
                if i < j {
                    if distance_2 <= roche_radius.powi(2) {
                        println!("\n");
                        panic!(
                            "[PANIC {} UTC] Particle {} was destroyed by particle {} due to close encounter!",
                            OffsetDateTime::now_utc()
                                .format(
                                    &format_description::parse(
                                        "[year].[month].[day] [hour]:[minute]:[second]"
                                    )
                                    .unwrap()
                                )
                                .unwrap(),
                            i,
                            j
                        );
                    }
                    // Check if particles are overlapping
                    if distance_2 <= (particle_a.radius + particle_b.radius).powi(2) {
                        println!("\n");
                        panic!(
                            "[PANIC {} UTC] Collision between particle {} and {}!",
                            OffsetDateTime::now_utc()
                                .format(
                                    &format_description::parse(
                                        "[year].[month].[day] [hour]:[minute]:[second]"
                                    )
                                    .unwrap()
                                )
                                .unwrap(),
                            i,
                            j
                        );
                    }
                    if MAX_DISTANCE_2 > 0.
                        && i == self.hosts.index.most_massive
                        && distance_2 > MAX_DISTANCE_2
                    {
                        println!("\n");
                        panic!(
                            "[PANIC {} UTC] Particle {} has been ejected!",
                            OffsetDateTime::now_utc()
                                .format(
                                    &format_description::parse(
                                        "[year].[month].[day] [hour]:[minute]:[second]"
                                    )
                                    .unwrap()
                                )
                                .unwrap(),
                            j
                        );
                    }
                }
                //////////////////////////////////////////////////////////////////////

                if ignore_terms == IgnoreGravityTerms::WHFastOne
                    && ((i == self.hosts.index.most_massive
                        && self.hosts.index.most_massive == 0
                        && j == 1)
                        || (i == self.hosts.index.most_massive
                            && self.hosts.index.most_massive > 0
                            && j == 0)
                        || (j == self.hosts.index.most_massive
                            && self.hosts.index.most_massive == 0
                            && i == 1)
                        || (j == self.hosts.index.most_massive
                            && self.hosts.index.most_massive > 0
                            && i == 0))
                {
                    // For WHFast Jacobi coordinates,
                    // ignore interaction between central body and first planet
                    continue;
                }
                if ignore_terms == IgnoreGravityTerms::WHFastTwo
                    && (i == self.hosts.index.most_massive || j == self.hosts.index.most_massive)
                {
                    // For WHFast democratic-heliocentric & WHDS coordinates
                    // completely ignore central body
                    continue;
                }

                // gravitational force equation for vector quantities
                // Fx = - G * ((m1 * m2) / r^3) * (x2 - x1)
                // Fy = - G * ((m1 * m2) / r^3) * (y2 - y1)
                // Fz = - G * ((m1 * m2) / r^3) * (z2 - z1)
                //
                // acceleration:
                // ax = Fx / m1 = - G * (m2 / r^3) * (x2 - x1)
                //let dx = particle_a.inertial_position.x - particle_b.inertial_position.x;
                //let dy = particle_a.inertial_position.y - particle_b.inertial_position.y;
                //let dz = particle_a.inertial_position.z - particle_b.inertial_position.z;
                //let r = sqrt!(dx*dx + dy*dy + dz*dz);
                //let r = sqrt!(distance_2 + softening);
                let distance = sqrt!(distance_2);
                let prefact = -G / (distance * distance * distance) * particle_b.mass;

                newtonian_inertial_acceleration.x += prefact * dx;
                newtonian_inertial_acceleration.y += prefact * dy;
                newtonian_inertial_acceleration.z += prefact * dz;

                ////// Compensated summation to improve the accuracy of additions that involve
                ////// one small and one large floating point number (already considered by IAS15)
                //////      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
                //////      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
                //// x
                //let acceleration_to_be_added = prefact * dx;
                //let corrected_acceleration_to_be_added = acceleration_to_be_added - newtonian_inertial_acceleration_error.x;
                //let added_acceleration = newtonian_inertial_acceleration.x + corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration_error.x = (added_acceleration - newtonian_inertial_acceleration.x) - corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration.x = added_acceleration;
                //// y
                //let acceleration_to_be_added = prefact * dy;
                //let corrected_acceleration_to_be_added = acceleration_to_be_added - newtonian_inertial_acceleration_error.y;
                //let added_acceleration = newtonian_inertial_acceleration.y + corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration_error.y = (added_acceleration - newtonian_inertial_acceleration.y) - corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration.y = added_acceleration;
                //// z
                //let acceleration_to_be_added = prefact * dz;
                //let corrected_acceleration_to_be_added = acceleration_to_be_added - newtonian_inertial_acceleration_error.z;
                //let added_acceleration = newtonian_inertial_acceleration.z + corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration_error.z = (added_acceleration - newtonian_inertial_acceleration.z) - corrected_acceleration_to_be_added;
                //newtonian_inertial_acceleration.z = added_acceleration;
            }
        }
        for (particle, newtonian_inertial_acceleration) in particles
            .iter_mut()
            .zip(newtonian_inertial_accelerations.iter())
        {
            particle.inertial_acceleration = *newtonian_inertial_acceleration;
        }
    }

    pub fn calculate_spin_and_evolving_quantities(&mut self, current_time: f64, evolution: bool) {
        let (particles, _) = self.particles.split_at_mut(self.n_particles);
        if evolution && self.consider_effects.evolution {
            // Evolve quantities (among others, it computes a new moment of inertia)
            evolution::calculate_particles_non_spin_dependent_evolving_quantities(
                current_time,
                particles,
                &mut self.particles_evolvers,
            );
        }
        common::calculate_spin(particles); // uses moment of inertia to compute spin
        if evolution && self.consider_effects.evolution {
            // Evolve quantities that require a computed spin
            evolution::calculate_particles_spin_dependent_evolving_quantities(
                current_time,
                particles,
                &mut self.particles_evolvers,
            );
        }
    }

    pub fn inertial_to_heliocentric(&mut self) {
        let (particles_left, particles_right) =
            self.particles.split_at_mut(self.hosts.index.most_massive);
        if let Some((host_particle, particles_right)) = particles_right.split_first_mut() {
            for particle in particles_left.iter_mut().chain(particles_right.iter_mut()) {
                particle.heliocentric_position.x =
                    particle.inertial_position.x - host_particle.inertial_position.x;
                particle.heliocentric_position.y =
                    particle.inertial_position.y - host_particle.inertial_position.y;
                particle.heliocentric_position.z =
                    particle.inertial_position.z - host_particle.inertial_position.z;
                particle.heliocentric_velocity.x =
                    particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.heliocentric_velocity.y =
                    particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.heliocentric_velocity.z =
                    particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.heliocentric_distance = sqrt!(
                    particle.heliocentric_position.x.powi(2)
                        + particle.heliocentric_position.y.powi(2)
                        + particle.heliocentric_position.z.powi(2)
                );
                particle.heliocentric_radial_velocity = (particle.heliocentric_position.x
                    * particle.heliocentric_velocity.x
                    + particle.heliocentric_position.y * particle.heliocentric_velocity.y
                    + particle.heliocentric_position.z * particle.heliocentric_velocity.z)
                    / particle.heliocentric_distance;
                particle.heliocentric_norm_velocity_vector_2 = (particle.heliocentric_velocity.x
                    - host_particle.heliocentric_velocity.x)
                    .powi(2)
                    + (particle.heliocentric_velocity.y - host_particle.heliocentric_velocity.y)
                        .powi(2)
                    + (particle.heliocentric_velocity.z - host_particle.heliocentric_velocity.z)
                        .powi(2);
                particle.heliocentric_norm_velocity_vector =
                    sqrt!(particle.heliocentric_norm_velocity_vector_2);
            }
            host_particle.heliocentric_position.zero();
            host_particle.heliocentric_velocity.zero();
            host_particle.heliocentric_distance = 0.;
            host_particle.heliocentric_radial_velocity = 0.;
            host_particle.heliocentric_norm_velocity_vector_2 = 0.;
            host_particle.heliocentric_norm_velocity_vector = 0.;
        }
    }

    pub fn initialize(&mut self, dangular_momentum_dt: bool, accelerations: bool) {
        let initialize_tides =
            self.consider_effects.tides && (accelerations || dangular_momentum_dt);
        let initialize_rotational_flattening =
            self.consider_effects.rotational_flattening && (dangular_momentum_dt || accelerations);
        let initialize_wind = dangular_momentum_dt && self.consider_effects.wind;
        let initialize_general_relativity = self.consider_effects.general_relativity
            && (accelerations
                || (dangular_momentum_dt
                    && self.general_relativity_implementation
                        == GeneralRelativityImplementation::Kidder1995));
        let initialize_disk = accelerations && self.consider_effects.disk;
        //
        let initialize = initialize_tides
            || initialize_rotational_flattening
            || initialize_wind
            || initialize_general_relativity
            || initialize_disk;
        if initialize {
            let (particles, _) = self.particles.split_at_mut(self.n_particles);
            //
            if initialize_wind {
                wind::initialize(particles);
            }
            //
            let (particles_left, particles_right) =
                particles.split_at_mut(self.hosts.index.most_massive);
            if let Some((host_particle, particles_right)) = particles_right.split_first_mut() {
                if initialize_tides && self.hosts.most_massive.tides {
                    tides::copy_heliocentric_coordinates(
                        host_particle,
                        particles_left,
                        particles_right,
                    );
                    tides::initialize(host_particle, particles_left, particles_right);
                }
                if initialize_rotational_flattening && self.hosts.most_massive.rotational_flattening
                {
                    rotational_flattening::copy_heliocentric_coordinates(
                        host_particle,
                        particles_left,
                        particles_right,
                    );
                    rotational_flattening::initialize(
                        host_particle,
                        particles_left,
                        particles_right,
                    );
                }
                if initialize_general_relativity && self.hosts.most_massive.general_relativity {
                    general_relativity::copy_heliocentric_coordinates(
                        host_particle,
                        particles_left,
                        particles_right,
                    );
                    general_relativity::initialize(host_particle, particles_left, particles_right);
                }
                if initialize_disk && self.hosts.most_massive.disk {
                    disk::copy_heliocentric_coordinates(
                        host_particle,
                        particles_left,
                        particles_right,
                    );
                    disk::initialize(host_particle, particles_left, particles_right);
                }
            }
            //
            if !self.hosts.most_massive.all {
                if initialize_tides && !self.hosts.most_massive.tides {
                    let (particles_left, particles_right) =
                        particles.split_at_mut(self.hosts.index.tides);
                    if let Some((tidal_host_particle, particles_right)) =
                        particles_right.split_first_mut()
                    {
                        tides::inertial_to_heliocentric_coordinates(
                            tidal_host_particle,
                            particles_left,
                            particles_right,
                        );
                        tides::initialize(tidal_host_particle, particles_left, particles_right);
                    }
                }
                if initialize_rotational_flattening
                    && !self.hosts.most_massive.rotational_flattening
                {
                    let (particles_left, particles_right) =
                        particles.split_at_mut(self.hosts.index.rotational_flattening);
                    if let Some((rotational_flattening_host_particle, particles_right)) =
                        particles_right.split_first_mut()
                    {
                        rotational_flattening::inertial_to_heliocentric_coordinates(
                            rotational_flattening_host_particle,
                            particles_left,
                            particles_right,
                        );
                        rotational_flattening::initialize(
                            rotational_flattening_host_particle,
                            particles_left,
                            particles_right,
                        );
                    }
                }
                if initialize_general_relativity && !self.hosts.most_massive.general_relativity {
                    let (particles_left, particles_right) =
                        particles.split_at_mut(self.hosts.index.general_relativity);
                    if let Some((general_relativity_host_particle, particles_right)) =
                        particles_right.split_first_mut()
                    {
                        general_relativity::inertial_to_heliocentric_coordinates(
                            general_relativity_host_particle,
                            particles_left,
                            particles_right,
                        );
                        general_relativity::initialize(
                            general_relativity_host_particle,
                            particles_left,
                            particles_right,
                        );
                    }
                }
                if initialize_disk && !self.hosts.most_massive.disk {
                    let (particles_left, particles_right) =
                        particles.split_at_mut(self.hosts.index.disk);
                    if let Some((disk_host_particle, particles_right)) =
                        particles_right.split_first_mut()
                    {
                        disk::inertial_to_heliocentric_coordinates(
                            disk_host_particle,
                            particles_left,
                            particles_right,
                        );
                        disk::initialize(disk_host_particle, particles_left, particles_right);
                    }
                }
            }
        }
        if accelerations {
            for particle in &mut self.particles[..self.n_particles] {
                particle.inertial_additional_acceleration.zero();
            }
        }
    }

    pub fn calculate_additional_effects(
        &mut self,
        current_time: f64,
        evolution: bool,
        dangular_momentum_dt: bool,
        accelerations: bool,
        ignored_gravity_terms: IgnoreGravityTerms,
    ) {
        self.initialize(dangular_momentum_dt, accelerations);
        self.calculate_spin_and_evolving_quantities(current_time, evolution); // Make sure we start with the good initial values

        let (particles, _) = self.particles.split_at_mut(self.n_particles);
        if dangular_momentum_dt && self.consider_effects.wind {
            wind::calculate_wind_factor(particles);
        }

        if self.consider_effects.tides || self.consider_effects.rotational_flattening {
            let (particles_left, particles_right) = particles.split_at_mut(self.hosts.index.tides);
            if let Some((tidal_host_particle, particles_right)) = particles_right.split_first_mut()
            {
                if (dangular_momentum_dt
                    && (self.consider_effects.tides || self.consider_effects.rotational_flattening))
                    || (accelerations
                        && (self.consider_effects.tides
                            || self.consider_effects.disk
                            || self.consider_effects.rotational_flattening
                            || self.consider_effects.general_relativity))
                {
                    {
                        //// calculate_orthogonal_components

                        tides::calculate_pair_dependent_scaled_dissipation_factors(
                            tidal_host_particle,
                            particles_left,
                            particles_right,
                            &mut self.pair_dependent_scaled_dissipation_factor,
                        ); // Needed by calculate_orthogonal_component_of_the_tidal_force and calculate_orthogonal_component_of_the_tidal_force if BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true)

                        if self.consider_effects.tides {
                            tides::calculate_orthogonal_component_of_the_tidal_force(
                                tidal_host_particle,
                                particles_left,
                                particles_right,
                                &mut self.pair_dependent_scaled_dissipation_factor,
                            );
                            tides::calculate_creep_coplanar_shapes(
                                tidal_host_particle,
                                particles_left,
                                particles_right,
                            );
                        }

                        if self.consider_effects.rotational_flattening {
                            rotational_flattening::calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(tidal_host_particle, particles_left, particles_right);
                            rotational_flattening::calculate_creep_coplanar_shapes(
                                tidal_host_particle,
                                particles_left,
                                particles_right,
                            );
                        }
                    }

                    if accelerations
                        && (self.consider_effects.tides
                            || self.consider_effects.rotational_flattening)
                    {
                        if self.consider_effects.tides {
                            tides::calculate_radial_component_of_the_tidal_force(
                                tidal_host_particle,
                                particles_left,
                                particles_right,
                                &mut self.pair_dependent_scaled_dissipation_factor,
                            ); // Needed for calculate_tidal_acceleration
                            tides::calculate_tidal_acceleration(
                                tidal_host_particle,
                                particles_left,
                                particles_right,
                            );
                        }

                        if self.consider_effects.rotational_flattening {
                            rotational_flattening::calculate_radial_component_of_the_force_induced_by_rotational_flattening(tidal_host_particle, particles_left, particles_right);
                            rotational_flattening::calculate_acceleration_induced_by_rotational_flattering(tidal_host_particle, particles_left, particles_right);
                        }
                    }

                    if dangular_momentum_dt
                        && (self.consider_effects.tides
                            || self.consider_effects.rotational_flattening)
                    {
                        // Not needed for additional accelerations
                        {
                            //// calculate_torques // Needed for dangular_momentum_dt
                            if self.consider_effects.tides {
                                tides::calculate_dangular_momentum_dt_due_to_tides(
                                    tidal_host_particle,
                                    particles_left,
                                    particles_right,
                                );
                            }

                            if self.consider_effects.rotational_flattening {
                                rotational_flattening::calculate_dangular_momentum_dt_induced_by_rotational_flattening(tidal_host_particle, particles_left, particles_right);
                            }
                        }
                    }
                }
            }
        }
        if accelerations && self.consider_effects.general_relativity {
            let (particles_left, particles_right) =
                particles.split_at_mut(self.hosts.index.general_relativity);
            if let Some((general_relativity_host_particle, particles_right)) =
                particles_right.split_first_mut()
            {
                if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) =
                    general_relativity_host_particle.general_relativity.effect
                {
                    match general_relativity_implementation {
                        GeneralRelativityImplementation::Kidder1995 => {
                            general_relativity::calculate_kidder1995_general_relativity_acceleration(general_relativity_host_particle, particles_left, particles_right);
                        }
                        GeneralRelativityImplementation::Anderson1975 => {
                            general_relativity::calculate_anderson1975_general_relativity_acceleration(general_relativity_host_particle, particles_left, particles_right, ignored_gravity_terms);
                        }
                        GeneralRelativityImplementation::Newhall1983 => {
                            general_relativity::calculate_newhall1983_general_relativity_acceleration(general_relativity_host_particle, particles_left, particles_right, ignored_gravity_terms);
                        }
                        GeneralRelativityImplementation::Disabled => {}
                    }
                }
            }
        }
        if accelerations && self.consider_effects.disk {
            let (particles_left, particles_right) = particles.split_at_mut(self.hosts.index.disk);
            if let Some((disk_host_particle, particles_right)) = particles_right.split_first_mut() {
                if let DiskEffect::CentralBody(_) = disk_host_particle.disk.effect {
                    disk::calculate_disk_interaction_acceleration(
                        current_time,
                        disk_host_particle,
                        particles_left,
                        particles_right,
                    );
                }
            }
        }
        if dangular_momentum_dt
            && (self.consider_effects.tides
                || self.consider_effects.rotational_flattening
                || (self.consider_effects.general_relativity
                    && self.general_relativity_implementation
                        == GeneralRelativityImplementation::Kidder1995)
                || self.consider_effects.wind)
        {
            self.calculate_dangular_momentum_dt();
        }

        if accelerations {
            self.add_additional_acceleration_corrections();
            //self.apply_acceleration_corrections();
            //println!("{:?}", self.particles[1].tides.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].rotational_flattening.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].general_relativity.parameters.output.acceleration);
            //println!("{:?}", self.particles[1].disk.parameters.output.acceleration);
            //panic!("*** {:?}", self.particles[1].inertial_acceleration);
        }
    }

    fn add_additional_acceleration_corrections(&mut self) {
        // Add the tidal+flattening+general relativity accelerations to the gravitational one (already computed)
        for particle in &mut self.particles[..self.n_particles] {
            if self.consider_effects.tides {
                particle.inertial_additional_acceleration.x +=
                    particle.tides.parameters.output.acceleration.x;
                particle.inertial_additional_acceleration.y +=
                    particle.tides.parameters.output.acceleration.y;
                particle.inertial_additional_acceleration.z +=
                    particle.tides.parameters.output.acceleration.z;
            }

            if self.consider_effects.disk {
                particle.inertial_additional_acceleration.x +=
                    particle.disk.parameters.output.acceleration.x;
                particle.inertial_additional_acceleration.y +=
                    particle.disk.parameters.output.acceleration.y;
                particle.inertial_additional_acceleration.z +=
                    particle.disk.parameters.output.acceleration.z;
            }

            if self.consider_effects.rotational_flattening {
                particle.inertial_additional_acceleration.x += particle
                    .rotational_flattening
                    .parameters
                    .output
                    .acceleration
                    .x;
                particle.inertial_additional_acceleration.y += particle
                    .rotational_flattening
                    .parameters
                    .output
                    .acceleration
                    .y;
                particle.inertial_additional_acceleration.z += particle
                    .rotational_flattening
                    .parameters
                    .output
                    .acceleration
                    .z;
            }

            if self.consider_effects.general_relativity {
                particle.inertial_additional_acceleration.x +=
                    particle.general_relativity.parameters.output.acceleration.x;
                particle.inertial_additional_acceleration.y +=
                    particle.general_relativity.parameters.output.acceleration.y;
                particle.inertial_additional_acceleration.z +=
                    particle.general_relativity.parameters.output.acceleration.z;
            }
        }
    }

    pub fn apply_acceleration_corrections(&mut self) {
        for particle in &mut self.particles[..self.n_particles] {
            particle.inertial_acceleration.x += particle.inertial_additional_acceleration.x;
            particle.inertial_acceleration.y += particle.inertial_additional_acceleration.y;
            particle.inertial_acceleration.z += particle.inertial_additional_acceleration.z;
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // TIDES
    fn calculate_dangular_momentum_dt(&mut self) {
        //// For compensated summation:
        //let mut dangular_momentum_dt_errors = [Axes{x:0., y:0., z:0. }; MAX_PARTICLES];
        //let (dangular_momentum_dt_errors, _) = dangular_momentum_dt_errors.split_at_mut(self.n_particles);

        //// For compensated summation:
        //for (particle, dangular_momentum_dt_error) in self.particles[..self.n_particles].iter_mut().zip(dangular_momentum_dt_errors.iter_mut()) {
        for particle in &mut self.particles[..self.n_particles] {
            // - Equation 25 from Bolmont et al. 2015
            particle.dangular_momentum_dt.x =
                particle.tides.parameters.output.dangular_momentum_dt.x
                    + particle
                        .rotational_flattening
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .x
                    + particle
                        .general_relativity
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .x
                    + particle.wind.parameters.output.dangular_momentum_dt.x;
            particle.dangular_momentum_dt.y =
                particle.tides.parameters.output.dangular_momentum_dt.y
                    + particle
                        .rotational_flattening
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .y
                    + particle
                        .general_relativity
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .y
                    + particle.wind.parameters.output.dangular_momentum_dt.y;
            particle.dangular_momentum_dt.z =
                particle.tides.parameters.output.dangular_momentum_dt.z
                    + particle
                        .rotational_flattening
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .z
                    + particle
                        .general_relativity
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .z
                    + particle.wind.parameters.output.dangular_momentum_dt.z;

            //// Compensated summation to improve the accuracy of additions that involve
            //// one small and one large floating point number (already considered by IAS15)
            ////      As cited by IAS15 paper: Kahan 1965; Higham 2002; Hairer et al. 2006
            ////      https://en.wikipedia.org/wiki/Kahan_summation_algorithm
            //for &dangular_momentum_dt_to_be_added in [particle.tides.parameters.output.dangular_momentum_dt, particle.rotational_flattening.parameters.output.dangular_momentum_dt, particle.general_relativity.parameters.output.dangular_momentum_dt, particle.wind.parameters.output.dangular_momentum_dt].iter() {
            //let corrected_dangular_momentum_dt_to_be_added = dangular_momentum_dt_to_be_added.x - dangular_momentum_dt_error.x;
            //let added_dangular_momentum_dt = particle.dangular_momentum_dt.x + corrected_dangular_momentum_dt_to_be_added;
            //dangular_momentum_dt_error.x = (added_dangular_momentum_dt - particle.dangular_momentum_dt.x) - corrected_dangular_momentum_dt_to_be_added;
            //particle.dangular_momentum_dt.x = added_dangular_momentum_dt;
            ////
            //let corrected_dangular_momentum_dt_to_be_added = dangular_momentum_dt_to_be_added.y - dangular_momentum_dt_error.y;
            //let added_dangular_momentum_dt = particle.dangular_momentum_dt.y + corrected_dangular_momentum_dt_to_be_added;
            //dangular_momentum_dt_error.y = (added_dangular_momentum_dt - particle.dangular_momentum_dt.y) - corrected_dangular_momentum_dt_to_be_added;
            //particle.dangular_momentum_dt.y = added_dangular_momentum_dt;
            ////
            //let corrected_dangular_momentum_dt_to_be_added = dangular_momentum_dt_to_be_added.z - dangular_momentum_dt_error.z;
            //let added_dangular_momentum_dt = particle.dangular_momentum_dt.z + corrected_dangular_momentum_dt_to_be_added;
            //dangular_momentum_dt_error.z = (added_dangular_momentum_dt - particle.dangular_momentum_dt.z) - corrected_dangular_momentum_dt_to_be_added;
            //particle.dangular_momentum_dt.z = added_dangular_momentum_dt;
            //}
        }
    }

    pub fn calculate_denergy_dt(&mut self) {
        let (particles_left, particles_right) =
            self.particles[..self.n_particles].split_at_mut(self.hosts.index.tides);
        if let Some((_, particles_right)) = particles_right.split_first_mut() {
            tides::calculate_denergy_dt(particles_left, particles_right);
        }
    }

    pub fn compute_total_energy(&self) -> f64 {
        let mut e_kin = 0.;
        let mut e_pot = 0.;
        let e_offset = 0.; // Energy offset due to collisions and ejections

        let (particles, _) = self.particles.split_at(self.n_particles);

        // Kinectic energy
        for particle in particles {
            e_kin += 0.5
                * particle.mass
                * (particle.heliocentric_velocity.x.powi(2)
                    + particle.heliocentric_velocity.y.powi(2)
                    + particle.heliocentric_velocity.z.powi(2));
        }
        // Gravitationl potential energy
        for (i, particle_a) in particles.iter().enumerate() {
            for particle_b in &particles[i + 1..] {
                let dx = particle_a.heliocentric_position.x - particle_b.heliocentric_position.x;
                let dy = particle_a.heliocentric_position.y - particle_b.heliocentric_position.y;
                let dz = particle_a.heliocentric_position.z - particle_b.heliocentric_position.z;
                e_pot -= particle_b.mass_g * particle_a.mass
                    / sqrt!(dx.powi(2) + dy.powi(2) + dz.powi(2));
            }
        }

        e_kin + e_pot + e_offset
    }

    pub fn compute_total_angular_momentum(&self) -> f64 {
        let mut total_angular_momentum = Axes::new(); // L
        for particle in &self.particles[..self.n_particles] {
            total_angular_momentum.x += particle.mass
                * (particle.heliocentric_position.y * particle.heliocentric_velocity.z
                    - particle.heliocentric_position.z * particle.heliocentric_velocity.y);
            total_angular_momentum.y += particle.mass
                * (particle.heliocentric_position.z * particle.heliocentric_velocity.x
                    - particle.heliocentric_position.x * particle.heliocentric_velocity.z);
            total_angular_momentum.z += particle.mass
                * (particle.heliocentric_position.x * particle.heliocentric_velocity.y
                    - particle.heliocentric_position.y * particle.heliocentric_velocity.x);
        }

        sqrt!(
            total_angular_momentum.x.powi(2)
                + total_angular_momentum.y.powi(2)
                + total_angular_momentum.z.powi(2)
        )
    }
}

fn get_center_of_mass_of_pair(
    center_of_mass_position: &mut Axes,
    center_of_mass_velocity: &mut Axes,
    center_of_mass_mass: f64,
    particle: &Particle,
) -> f64 {
    center_of_mass_position.x = center_of_mass_position.x * center_of_mass_mass
        + particle.heliocentric_position.x * particle.mass;
    center_of_mass_position.y = center_of_mass_position.y * center_of_mass_mass
        + particle.heliocentric_position.y * particle.mass;
    center_of_mass_position.z = center_of_mass_position.z * center_of_mass_mass
        + particle.heliocentric_position.z * particle.mass;
    center_of_mass_velocity.x = center_of_mass_velocity.x * center_of_mass_mass
        + particle.heliocentric_velocity.x * particle.mass;
    center_of_mass_velocity.y = center_of_mass_velocity.y * center_of_mass_mass
        + particle.heliocentric_velocity.y * particle.mass;
    center_of_mass_velocity.z = center_of_mass_velocity.z * center_of_mass_mass
        + particle.heliocentric_velocity.z * particle.mass;

    let new_center_of_mass_mass = center_of_mass_mass + particle.mass;
    if new_center_of_mass_mass > 0. {
        center_of_mass_position.x /= new_center_of_mass_mass;
        center_of_mass_position.y /= new_center_of_mass_mass;
        center_of_mass_position.z /= new_center_of_mass_mass;
        center_of_mass_velocity.x /= new_center_of_mass_mass;
        center_of_mass_velocity.y /= new_center_of_mass_mass;
        center_of_mass_velocity.z /= new_center_of_mass_mass;
    }
    new_center_of_mass_mass
}

pub fn calculate_center_of_mass(particles: &[Particle]) -> (Axes, Axes) {
    // Compute center of mass
    let mut center_of_mass_position = Axes::new();
    let mut center_of_mass_velocity = Axes::new();
    let mut center_of_mass_mass = 0.;

    for particle in particles {
        center_of_mass_mass = get_center_of_mass_of_pair(
            &mut center_of_mass_position,
            &mut center_of_mass_velocity,
            center_of_mass_mass,
            particle,
        );
    }

    (center_of_mass_position, center_of_mass_velocity)
}

fn disable_unnecessary_effects(consider_effects: &mut ConsiderEffects, particles: &[Particle]) {
    let mut found_central_body_tides = false;
    let mut found_central_body_rotational_flattening = false;
    let mut found_central_body_general_relativity = false;
    let mut found_central_body_disk = false;
    let mut found_wind = false;
    let mut found_evolving_body = false;
    for particle in particles {
        if let TidesEffect::CentralBody(_) = particle.tides.effect {
            assert!(
                !found_central_body_tides,
                "[PANIC {} UTC] Only one central body is allowed for tidal effect!",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            found_central_body_tides = true;
        }
        if let RotationalFlatteningEffect::CentralBody(_) = particle.rotational_flattening.effect {
            assert!(
                !found_central_body_rotational_flattening,
                "[PANIC {} UTC] Only one central body is allowed for rotational flattening effects!",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            found_central_body_rotational_flattening = true;
        }
        if let GeneralRelativityEffect::CentralBody(implementation) =
            particle.general_relativity.effect
        {
            if implementation != GeneralRelativityImplementation::Disabled {
                assert!(
                    !found_central_body_general_relativity,
                    "[PANIC {} UTC] Only one central body is allowed for general relativity effects!",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap()
                );
                found_central_body_general_relativity = true;
            }
        }
        if let DiskEffect::CentralBody(_) = particle.disk.effect {
            assert!(
                !found_central_body_disk,
                "[PANIC {} UTC] Only one central body is allowed for disk effects!",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            found_central_body_disk = true;
        }
        if WindEffect::Disabled != particle.wind.effect {
            found_wind = true;
        }
        if EvolutionType::NonEvolving != particle.evolution {
            found_evolving_body = true;
        }
    }
    if consider_effects.tides && !found_central_body_tides {
        println!(
            "[WARNING {} UTC] Disabled tides because no central host was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.tides = false;
    }
    if consider_effects.rotational_flattening && !found_central_body_rotational_flattening {
        println!(
            "[WARNING {} UTC] Disabled rotational flattening because no central host was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.rotational_flattening = false;
    }
    if consider_effects.general_relativity && !found_central_body_general_relativity {
        println!(
            "[WARNING {} UTC] Disabled general relativity because no central host was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.general_relativity = false;
    }
    if consider_effects.disk && !found_central_body_disk {
        println!(
            "[WARNING {} UTC] Disabled disk because no central host was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.disk = false;
    }
    if consider_effects.wind && !found_wind {
        println!(
            "[WARNING {} UTC] Disabled wind because no wind was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.wind = false;
    }
    if consider_effects.evolution && !found_evolving_body {
        println!(
            "[WARNING {} UTC] Disabled evolution because no evolving body was included!",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        consider_effects.evolution = false;
    }
}

fn check_effects_vs_central_and_orbiting(
    particles: &[Particle],
    consider_effects: &ConsiderEffects,
) {
    let mut found_tides_central_body = false;
    let mut found_tides_orbiting_body = false;
    for (i, particle) in particles.iter().enumerate() {
        if let TidesEffect::CentralBody(_) = particle.tides.effect {
            found_tides_central_body = true;
            if !consider_effects.tides {
                println!(
                    "[WARNING {} UTC] Particle {} has tidal effect (central body) but the tidal effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
        if let TidesEffect::OrbitingBody(_) = particle.tides.effect {
            found_tides_orbiting_body = true;
            if !consider_effects.tides {
                println!(
                    "[WARNING {} UTC] Particle {} has tidal effect (orbiting body) but the tidal effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.tides {
        if !found_tides_central_body {
            println!(
                "[INFO {} UTC] No central body for tidal effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        if !found_tides_orbiting_body {
            println!(
                "[INFO {} UTC] No orbiting body for tidal effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
    }

    let mut found_rotational_flattening_central_body = false;
    let mut found_rotational_flattening_orbiting_body = false;
    for (i, particle) in particles.iter().enumerate() {
        if let RotationalFlatteningEffect::CentralBody(_) = particle.rotational_flattening.effect {
            found_rotational_flattening_central_body = true;
            if !consider_effects.rotational_flattening {
                println!(
                    "[WARNING {} UTC] Particle {} has rotational flattening effect (central body) but the rotational flattening effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
        if let RotationalFlatteningEffect::OrbitingBody(_) = particle.rotational_flattening.effect {
            found_rotational_flattening_orbiting_body = true;
            if !consider_effects.rotational_flattening {
                println!(
                    "[WARNING {} UTC] Particle {} has rotatial flattening effect (orbiting body) but the rotatial flattening effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.rotational_flattening {
        if !found_rotational_flattening_central_body {
            println!(
                "[INFO {} UTC] No central body for rotational flattening effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        if !found_rotational_flattening_orbiting_body {
            println!(
                "[INFO {} UTC] No orbiting body for rotational flattening effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
    }

    let mut found_general_relativity_central_body = false;
    let mut found_general_relativity_orbiting_body = false;
    for (i, particle) in particles.iter().enumerate() {
        if let GeneralRelativityEffect::CentralBody(implementation) =
            particle.general_relativity.effect
        {
            if implementation != GeneralRelativityImplementation::Disabled {
                found_general_relativity_central_body = true;
                if !consider_effects.general_relativity {
                    println!(
                        "[WARNING {} UTC] Particle {} has general relativity effect (central body) but the general relativity effect is disabled for this simulation",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap(),
                        i
                    );
                }
            }
        }
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            found_general_relativity_orbiting_body = true;
            if !consider_effects.general_relativity {
                println!(
                    "[WARNING {} UTC] Particle {} has general relativity effect (orbiting body) but the general relativity effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.general_relativity {
        if !found_general_relativity_central_body {
            println!(
                "[INFO {} UTC] No central body for general relativity effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        if !found_general_relativity_orbiting_body {
            println!(
                "[INFO {} UTC] No orbiting body for general relativity effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
    }

    let mut found_disk_central_body = false;
    let mut found_disk_orbiting_body = false;
    for (i, particle) in particles.iter().enumerate() {
        if let DiskEffect::CentralBody(_disk) = particle.disk.effect {
            found_disk_central_body = true;
            if !consider_effects.disk {
                println!(
                    "[WARNING {} UTC] Particle {} has disk effect (central body) but the disk effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
        if let DiskEffect::OrbitingBody = particle.disk.effect {
            found_disk_orbiting_body = true;
            if !consider_effects.disk {
                println!(
                    "[WARNING {} UTC] Particle {} has disk effect (orbiting body) but the disk effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.disk {
        if !found_disk_central_body {
            println!(
                "[INFO {} UTC] No central body for disk effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        if !found_disk_orbiting_body {
            println!(
                "[INFO {} UTC] No orbiting body for disk effects",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
    }

    let mut found_wind = false;
    for (i, particle) in particles.iter().enumerate() {
        if let WindEffect::Interaction = particle.wind.effect {
            found_wind = true;
            if !consider_effects.wind {
                println!(
                    "[WARNING {} UTC] Particle {} has wind effect (central body) but the wind effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.wind && !found_wind {
        println!(
            "[INFO {} UTC] No wind effects",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
    }

    let mut found_evolution = false;
    for (i, particle) in particles.iter().enumerate() {
        if particle.evolution != EvolutionType::NonEvolving {
            found_evolution = true;
            if !consider_effects.evolution {
                println!(
                    "[WARNING {} UTC] Particle {} has evolution effect but the evolution effect is disabled for this simulation",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    i
                );
            }
        }
    }
    if consider_effects.evolution && !found_evolution {
        println!(
            "[INFO {} UTC] No evolution effects",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
    }
}

fn find_indices(particles: &[Particle], consider_effects: &ConsiderEffects) -> Hosts {
    // Most massive particle
    let mut most_massive_particle_index = MAX_PARTICLES + 1;
    let mut max_mass_found = 0.;
    for (i, particle) in particles.iter().enumerate() {
        if particle.mass > max_mass_found {
            max_mass_found = particle.mass;
            most_massive_particle_index = i;
        }
    }

    // Particle that is the main one for general relativity effects and WHFast symplectic integration
    let mut general_relativity_host_particle_index = MAX_PARTICLES + 1;
    if consider_effects.general_relativity {
        for (i, particle) in particles.iter().enumerate() {
            if let GeneralRelativityEffect::CentralBody(implementation) =
                particle.general_relativity.effect
            {
                if implementation != GeneralRelativityImplementation::Disabled {
                    if general_relativity_host_particle_index == MAX_PARTICLES + 1 {
                        general_relativity_host_particle_index = i;
                    } else {
                        panic!("Only one central body is allowed for general relativity effects!");
                    }
                }
            }
        }
        assert!(
            (most_massive_particle_index == general_relativity_host_particle_index),
            "The central body for General Relativity should be the most massive one!"
        );
    }

    // Particle that is the main one for tidal effects
    let mut tidal_host_particle_index = MAX_PARTICLES + 1;
    if consider_effects.tides {
        for (i, particle) in particles.iter().enumerate() {
            if let TidesEffect::CentralBody(_) = particle.tides.effect {
                if tidal_host_particle_index == MAX_PARTICLES + 1 {
                    tidal_host_particle_index = i;
                } else {
                    panic!("Only one central body is allowed for tidal effects!");
                }
            }
        }
    }

    // Particle that is the main one for rotational flatting effects
    let mut rotational_flattening_host_particle_index = MAX_PARTICLES + 1;
    if consider_effects.rotational_flattening {
        for (i, particle) in particles.iter().enumerate() {
            if let RotationalFlatteningEffect::CentralBody(_) =
                particle.rotational_flattening.effect
            {
                if rotational_flattening_host_particle_index == MAX_PARTICLES + 1 {
                    rotational_flattening_host_particle_index = i;
                } else {
                    panic!("Only one central body is allowed for rotational flattening effects!");
                }
            }
        }
    }

    if consider_effects.tides
        && consider_effects.rotational_flattening
        && tidal_host_particle_index != rotational_flattening_host_particle_index
    {
        panic!("The central body for tidal & rotational flattening effects needs to be the same!");
    } else if !consider_effects.tides && consider_effects.rotational_flattening {
        tidal_host_particle_index = rotational_flattening_host_particle_index;
    }

    // Find the position of the particle that has the disk if any
    let mut disk_host_particle_index = MAX_PARTICLES + 1;
    if consider_effects.disk {
        for (i, particle) in particles.iter().enumerate() {
            if let DiskEffect::CentralBody(disk) = particle.disk.effect {
                if disk.inner_edge_distance != 0. || disk.outer_edge_distance != 0. {
                    if disk_host_particle_index == MAX_PARTICLES + 1 {
                        disk_host_particle_index = i;
                    } else {
                        panic!("Only one central body with a disk is allowed!");
                    }
                }
            }
        }
    }

    let tidal_host_particle_is_the_most_massive =
        consider_effects.tides && most_massive_particle_index == tidal_host_particle_index;
    let rotational_flattening_host_particle_is_the_most_massive = consider_effects
        .rotational_flattening
        && most_massive_particle_index == rotational_flattening_host_particle_index;
    let general_relativity_host_particle_is_the_most_massive = consider_effects.general_relativity
        && most_massive_particle_index == general_relativity_host_particle_index;
    let disk_host_particle_is_the_most_massive =
        consider_effects.disk && most_massive_particle_index == disk_host_particle_index;
    let all_same_particle_indices = (rotational_flattening_host_particle_is_the_most_massive
        || !consider_effects.rotational_flattening)
        && (tidal_host_particle_is_the_most_massive || !consider_effects.tides);

    Hosts {
        index: HostIndices {
            most_massive: most_massive_particle_index,
            tides: tidal_host_particle_index,
            rotational_flattening: rotational_flattening_host_particle_index,
            general_relativity: general_relativity_host_particle_index,
            disk: disk_host_particle_index,
        },
        most_massive: HostMostMassive {
            all: all_same_particle_indices,
            tides: tidal_host_particle_is_the_most_massive,
            rotational_flattening: rotational_flattening_host_particle_is_the_most_massive,
            general_relativity: general_relativity_host_particle_is_the_most_massive,
            disk: disk_host_particle_is_the_most_massive,
        },
    }
}
