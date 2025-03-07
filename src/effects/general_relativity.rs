use crate::constants::{DBL_EPSILON_2, G, MAX_PARTICLES, SPEED_OF_LIGHT_2};
use crate::particles::universe::IgnoreGravityTerms;
use crate::{Axes, Particle};
use itertools::izip;
use serde::{Deserialize, Serialize};
use std::iter;
use time;
use time::{OffsetDateTime, format_description};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct GeneralRelativityParticleInternalParameters {
    pub distance: f64,
    pub radial_velocity: f64,
    pub norm_velocity_vector: f64,
    pub norm_velocity_vector_2: f64,
    pub factor: f64,
}

impl GeneralRelativityParticleInternalParameters {
    pub fn zero(&mut self) {
        self.distance = 0.;
        self.radial_velocity = 0.;
        self.norm_velocity_vector = 0.;
        self.norm_velocity_vector_2 = 0.;
        self.factor = 0.;
    }
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleParameters {
    pub internal: GeneralRelativityParticleInternalParameters,
    pub output: GeneralRelativityParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativityParticleCoordinates {
    // Positions/velocities in a heliocentric frame
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum GeneralRelativityImplementation {
    Kidder1995,   // MercuryT
    Anderson1975, // REBOUNDx gr
    Newhall1983,  // REBOUNDx gr full
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum GeneralRelativityEffect {
    CentralBody(GeneralRelativityImplementation),
    OrbitingBody,
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct GeneralRelativity {
    pub effect: GeneralRelativityEffect,
    pub parameters: GeneralRelativityParticleParameters,
    pub coordinates: GeneralRelativityParticleCoordinates,
}

impl GeneralRelativity {
    pub fn new(effect: GeneralRelativityEffect) -> GeneralRelativity {
        GeneralRelativity {
            effect,
            parameters: GeneralRelativityParticleParameters {
                internal: GeneralRelativityParticleInternalParameters::default(),
                output: GeneralRelativityParticleOutputParameters {
                    acceleration: Axes::new(),
                    dangular_momentum_dt: Axes::new(),
                },
            },
            coordinates: GeneralRelativityParticleCoordinates {
                position: Axes::new(),
                velocity: Axes::new(),
            },
        }
    }
}

pub fn initialize(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) =
        host_particle.general_relativity.effect
    {
        if general_relativity_implementation != GeneralRelativityImplementation::Newhall1983
            && general_relativity_implementation != GeneralRelativityImplementation::Disabled
        {
            host_particle.general_relativity.parameters.internal.factor = 0.;
            host_particle
                .general_relativity
                .parameters
                .output
                .acceleration
                .zero();
            host_particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .zero();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.parameters.internal.factor = host_particle.mass_g
                        * particle.mass_g
                        / (host_particle.mass_g + particle.mass_g).powi(2);
                    particle
                        .general_relativity
                        .parameters
                        .output
                        .acceleration
                        .zero();
                    particle
                        .general_relativity
                        .parameters
                        .output
                        .dangular_momentum_dt
                        .zero();
                }
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) =
        host_particle.general_relativity.effect
    {
        if general_relativity_implementation != GeneralRelativityImplementation::Disabled {
            // Inertial to Heliocentric positions/velocities
            host_particle.general_relativity.coordinates.position.zero();
            host_particle.general_relativity.coordinates.velocity.zero();
            // Clears all except for internal.factor
            host_particle.general_relativity.parameters.internal.zero();
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.coordinates.position = particle.inertial_position;
                    particle
                        .general_relativity
                        .coordinates
                        .position
                        .sub(&host_particle.inertial_position);

                    particle.general_relativity.coordinates.velocity = particle.inertial_velocity;
                    particle
                        .general_relativity
                        .coordinates
                        .velocity
                        .sub(&host_particle.inertial_velocity);

                    particle.general_relativity.parameters.internal.distance =
                        particle.general_relativity.coordinates.position.norm();

                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .radial_velocity = particle
                        .general_relativity
                        .coordinates
                        .position
                        .dot(&particle.general_relativity.coordinates.velocity)
                        / particle.general_relativity.parameters.internal.distance;

                    let mut tmp = particle.general_relativity.coordinates.velocity;
                    tmp.sub(&host_particle.general_relativity.coordinates.velocity);
                    let norm = tmp.norm();
                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .norm_velocity_vector = norm;

                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .norm_velocity_vector_2 = norm.powi(2);
                }
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let GeneralRelativityEffect::CentralBody(general_relativity_implementation) =
        host_particle.general_relativity.effect
    {
        if general_relativity_implementation != GeneralRelativityImplementation::Disabled {
            host_particle.general_relativity.coordinates.position =
                host_particle.heliocentric_position;
            host_particle.general_relativity.coordinates.velocity =
                host_particle.heliocentric_velocity;
            host_particle
                .general_relativity
                .parameters
                .internal
                .distance = host_particle.heliocentric_distance;
            host_particle
                .general_relativity
                .parameters
                .internal
                .radial_velocity = host_particle.heliocentric_radial_velocity;
            host_particle
                .general_relativity
                .parameters
                .internal
                .norm_velocity_vector = host_particle.heliocentric_norm_velocity_vector;
            host_particle
                .general_relativity
                .parameters
                .internal
                .norm_velocity_vector_2 = host_particle.heliocentric_norm_velocity_vector_2;
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
                    particle.general_relativity.coordinates.position =
                        particle.heliocentric_position;
                    particle.general_relativity.coordinates.velocity =
                        particle.heliocentric_velocity;
                    particle.general_relativity.parameters.internal.distance =
                        particle.heliocentric_distance;
                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .radial_velocity = particle.heliocentric_radial_velocity;
                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .norm_velocity_vector = particle.heliocentric_norm_velocity_vector;
                    particle
                        .general_relativity
                        .parameters
                        .internal
                        .norm_velocity_vector_2 = particle.heliocentric_norm_velocity_vector_2;
                }
            }
        }
    }
}

pub fn calculate_kidder1995_general_relativity_acceleration(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    calculate_kidder1995_first_order_general_relativity_acceleration(
        host_particle,
        particles,
        more_particles,
    );
    calculate_kidder1995_second_order_general_relativity_acceleration(
        host_particle,
        particles,
        more_particles,
    );
    calculate_kidder1995_spin_orbit_general_relativity_acceleration_and_dangular_momentum_dt(
        host_particle,
        particles,
        more_particles,
    );
}

fn calculate_kidder1995_first_order_general_relativity_acceleration(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    let mut sum_total_general_relativity_acceleration = Axes::new();

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            // Radial part of the GR force (Kidder 1995, Mardling & Lin 2002)
            // - Equation 11 from Bolmont et al. 2015
            let star_planet_mass_g = host_particle.mass_g + particle.mass_g;
            let distance_2 = particle
                .general_relativity
                .parameters
                .internal
                .distance
                .powi(2);
            let radial_velocity_2 = particle
                .general_relativity
                .parameters
                .internal
                .radial_velocity
                .powi(2);
            let radial_component_of_the_general_relativity_force = -star_planet_mass_g
                / (distance_2 * SPEED_OF_LIGHT_2)
                * ((1.0 + 3.0 * particle.general_relativity.parameters.internal.factor)
                    * particle
                        .general_relativity
                        .parameters
                        .internal
                        .norm_velocity_vector_2
                    - 2.0
                        * (2.0 + particle.general_relativity.parameters.internal.factor)
                        * star_planet_mass_g
                        / particle.general_relativity.parameters.internal.distance
                    - 1.5
                        * particle.general_relativity.parameters.internal.factor
                        * radial_velocity_2);
            //println!("Radial component GR force {:e}", radial_component_of_the_general_relativity_force);
            // Orthoradial part of the GR force
            // - Equation 11 from Bolmont et al. 2015
            let orthogonal_component_of_the_general_relativity_force = star_planet_mass_g
                / (distance_2 * SPEED_OF_LIGHT_2)
                * 2.0
                * (2.0 - particle.general_relativity.parameters.internal.factor)
                * particle
                    .general_relativity
                    .parameters
                    .internal
                    .radial_velocity
                * particle
                    .general_relativity
                    .parameters
                    .internal
                    .norm_velocity_vector;
            //println!("Ortho component GR force {:e}", orthogonal_component_of_the_general_relativity_force);
            // Total General Relativity force
            // - Equation 10 from Bolmont et al. 2015
            let mut total_general_relativity_acceleration = Axes::from(
                radial_component_of_the_general_relativity_force
                    * particle.general_relativity.coordinates.position.x
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force
                        * particle.general_relativity.coordinates.velocity.x
                        / particle
                            .general_relativity
                            .parameters
                            .internal
                            .norm_velocity_vector,
                radial_component_of_the_general_relativity_force
                    * particle.general_relativity.coordinates.position.y
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force
                        * particle.general_relativity.coordinates.velocity.y
                        / particle
                            .general_relativity
                            .parameters
                            .internal
                            .norm_velocity_vector,
                radial_component_of_the_general_relativity_force
                    * particle.general_relativity.coordinates.position.z
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_general_relativity_force
                        * particle.general_relativity.coordinates.velocity.z
                        / particle
                            .general_relativity
                            .parameters
                            .internal
                            .norm_velocity_vector,
            );

            // - Equation 19 from Bolmont et al. 2015 (first term)
            particle.general_relativity.parameters.output.acceleration =
                total_general_relativity_acceleration;

            total_general_relativity_acceleration.mul(particle.mass / host_particle.mass);
            sum_total_general_relativity_acceleration.add(&total_general_relativity_acceleration);

            //println!("GR force {:e} {:e} {:e}", total_general_relativity_acceleration.x,
            //total_general_relativity_acceleration.y, total_general_relativity_acceleration.z);
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
    //particle.general_relativity.parameters.output.acceleration.x += sum_total_general_relativity_acceleration.x;
    //particle.general_relativity.parameters.output.acceleration.y += sum_total_general_relativity_acceleration.y;
    //particle.general_relativity.parameters.output.acceleration.z += sum_total_general_relativity_acceleration.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    sum_total_general_relativity_acceleration.negate();
    host_particle
        .general_relativity
        .parameters
        .output
        .acceleration = sum_total_general_relativity_acceleration;
}

fn calculate_kidder1995_second_order_general_relativity_acceleration(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    // 2nd order Post-Newtonian
    let mut sum_total_second_order_general_relativity_acceleration = Axes::new();

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let star_planet_mass_g = host_particle.mass_g + particle.mass_g;
            let distance_2 = particle
                .general_relativity
                .parameters
                .internal
                .distance
                .powi(2);
            let norm_velocity_vector_2 = particle
                .general_relativity
                .parameters
                .internal
                .norm_velocity_vector_2;
            let norm_velocity_vector_4 = norm_velocity_vector_2.powi(2);
            let radial_velocity_2 = particle
                .general_relativity
                .parameters
                .internal
                .radial_velocity
                .powi(2);
            let radial_velocity_4 = radial_velocity_2.powi(2);
            let general_relativity_factor_2 = particle
                .general_relativity
                .parameters
                .internal
                .factor
                .powi(2);

            // Radial part of the GR force (Kidder 1995, equation 2.2d)
            let radial_component_of_the_second_order_general_relativity_acceleration =
                -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * (3.0 / 4.0
                        * (12.0 + 29.0 * particle.general_relativity.parameters.internal.factor)
                        * (star_planet_mass_g.powi(2) / distance_2)
                        + particle.general_relativity.parameters.internal.factor
                            * (3.0 - 4.0 * particle.general_relativity.parameters.internal.factor)
                            * norm_velocity_vector_4
                        + 15.0 / 8.0
                            * particle.general_relativity.parameters.internal.factor
                            * (1.0 - 3.0 * particle.general_relativity.parameters.internal.factor)
                            * radial_velocity_4
                        - 3.0 / 2.0
                            * particle.general_relativity.parameters.internal.factor
                            * (3.0 - 4.0 * particle.general_relativity.parameters.internal.factor)
                            * radial_velocity_2
                            * norm_velocity_vector_2
                        - 0.5
                            * particle.general_relativity.parameters.internal.factor
                            * (13.0
                                - 4.0 * particle.general_relativity.parameters.internal.factor)
                            * (star_planet_mass_g
                                / particle.general_relativity.parameters.internal.distance)
                            * norm_velocity_vector_2
                        - (2.0
                            + 25.0 * particle.general_relativity.parameters.internal.factor
                            + 2.0 * general_relativity_factor_2)
                            * (star_planet_mass_g
                                / particle.general_relativity.parameters.internal.distance)
                            * radial_velocity_2);

            let orthogonal_component_of_the_second_order_general_relativity_acceleration =
                -star_planet_mass_g / (distance_2 * SPEED_OF_LIGHT_2)
                    * (-0.5)
                    * particle
                        .general_relativity
                        .parameters
                        .internal
                        .radial_velocity
                    * (particle.general_relativity.parameters.internal.factor
                        * (15.0 + 4.0 * particle.general_relativity.parameters.internal.factor)
                        * norm_velocity_vector_2
                        - (4.0
                            + 41.0 * particle.general_relativity.parameters.internal.factor
                            + 8.0 * general_relativity_factor_2)
                            * (star_planet_mass_g
                                / particle.general_relativity.parameters.internal.distance)
                        - 3.0
                            * particle.general_relativity.parameters.internal.factor
                            * (3.0 + 2.0 * particle.general_relativity.parameters.internal.factor)
                            * radial_velocity_2);

            let mut total_second_order_general_relativity_acceleration = Axes::from(
                radial_component_of_the_second_order_general_relativity_acceleration
                    * particle.general_relativity.coordinates.position.x
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration
                        * particle.general_relativity.coordinates.velocity.x,
                radial_component_of_the_second_order_general_relativity_acceleration
                    * particle.general_relativity.coordinates.position.y
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration
                        * particle.general_relativity.coordinates.velocity.y,
                radial_component_of_the_second_order_general_relativity_acceleration
                    * particle.general_relativity.coordinates.position.z
                    / particle.general_relativity.parameters.internal.distance
                    + orthogonal_component_of_the_second_order_general_relativity_acceleration
                        * particle.general_relativity.coordinates.velocity.z,
            );

            particle
                .general_relativity
                .parameters
                .output
                .acceleration
                .add(&total_second_order_general_relativity_acceleration);

            total_second_order_general_relativity_acceleration
                .mul(particle.mass / host_particle.mass);
            sum_total_second_order_general_relativity_acceleration
                .add(&total_second_order_general_relativity_acceleration);
            //println!("a {} {} {}", total_second_order_general_relativity_acceleration.x, total_second_order_general_relativity_acceleration.y, total_second_order_general_relativity_acceleration.z);
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
    //particle.general_relativity.parameters.output.acceleration.x += sum_total_second_order_general_relativity_acceleration.x;
    //particle.general_relativity.parameters.output.acceleration.y += sum_total_second_order_general_relativity_acceleration.y;
    //particle.general_relativity.parameters.output.acceleration.z += sum_total_second_order_general_relativity_acceleration.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    sum_total_second_order_general_relativity_acceleration.negate();
    host_particle
        .general_relativity
        .parameters
        .output
        .acceleration
        .add(&sum_total_second_order_general_relativity_acceleration);
}

fn calculate_kidder1995_spin_orbit_general_relativity_acceleration_and_dangular_momentum_dt(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    // Equation 5 from https://arxiv.org/pdf/1102.5192.pdf
    // Spin effects are known for the dominant relativistic spin-orbit coupling term at 1.5PN
    // https://arxiv.org/pdf/gr-qc/0202016.pdf
    // Spin in Kidder is defined as angular_momentum
    let mut star_angular_momentum = host_particle.spin;
    star_angular_momentum.mul(host_particle.moment_of_inertia);
    host_particle
        .general_relativity
        .parameters
        .output
        .dangular_momentum_dt
        .zero();
    let mut sum_total_general_relativity_spin_orbit_acceleration = Axes::new();

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            // - Equation 2.2c from Kidder 1995
            let star_planet_mass = host_particle.mass + particle.mass;
            let star_planet_diff_mass = host_particle.mass - particle.mass;
            let mass_factor = star_planet_diff_mass / star_planet_mass;

            // Spin in Kidder is defined as angular_momentum
            let mut particle_angular_momentum = particle.spin;
            particle_angular_momentum.mul(particle.moment_of_inertia);

            let mut particle_normalized_position = particle.general_relativity.coordinates.position;
            particle_normalized_position
                .div(particle.general_relativity.parameters.internal.distance);

            let mass_spin_factor = Axes::from(
                mass_factor
                    * star_planet_mass
                    * (particle_angular_momentum.x / particle.mass
                        - star_angular_momentum.x / host_particle.mass),
                mass_factor
                    * star_planet_mass
                    * (particle_angular_momentum.y / particle.mass
                        - star_angular_momentum.y / host_particle.mass),
                mass_factor
                    * star_planet_mass
                    * (particle_angular_momentum.z / particle.mass
                        - star_angular_momentum.z / host_particle.mass),
            );

            let element1 = Axes::from(
                6. * particle_normalized_position.x
                    * ((particle_normalized_position.y
                        * particle.general_relativity.coordinates.velocity.z
                        - particle_normalized_position.z
                            * particle.general_relativity.coordinates.velocity.y)
                        * (2. * (star_angular_momentum.x + particle_angular_momentum.x)
                            + mass_spin_factor.x)),
                6. * particle_normalized_position.y
                    * ((particle_normalized_position.z
                        * particle.general_relativity.coordinates.velocity.x
                        - particle_normalized_position.x
                            * particle.general_relativity.coordinates.velocity.z)
                        * (2. * (star_angular_momentum.y + particle_angular_momentum.y)
                            + mass_spin_factor.y)),
                6. * particle_normalized_position.z
                    * ((particle_normalized_position.x
                        * particle.general_relativity.coordinates.velocity.y
                        - particle_normalized_position.y
                            * particle.general_relativity.coordinates.velocity.x)
                        * (2. * (star_angular_momentum.z + particle_angular_momentum.z)
                            + mass_spin_factor.z)),
            );

            let element7s = Axes::from(
                7. * (star_angular_momentum.x + particle_angular_momentum.x)
                    + 3. * mass_spin_factor.x,
                7. * (star_angular_momentum.y + particle_angular_momentum.y)
                    + 3. * mass_spin_factor.y,
                7. * (star_angular_momentum.z + particle_angular_momentum.z)
                    + 3. * mass_spin_factor.z,
            );
            let element2 = Axes::from(
                particle.general_relativity.coordinates.velocity.y * element7s.z
                    - particle.general_relativity.coordinates.velocity.z * element7s.y,
                particle.general_relativity.coordinates.velocity.z * element7s.x
                    - particle.general_relativity.coordinates.velocity.x * element7s.z,
                particle.general_relativity.coordinates.velocity.x * element7s.y
                    - particle.general_relativity.coordinates.velocity.y * element7s.x,
            );

            let element3s = Axes::from(
                3. * (star_angular_momentum.x + particle_angular_momentum.x) + mass_spin_factor.x,
                3. * (star_angular_momentum.y + particle_angular_momentum.y) + mass_spin_factor.y,
                3. * (star_angular_momentum.z + particle_angular_momentum.z) + mass_spin_factor.z,
            );
            let element3 = Axes::from(
                3. * particle
                    .general_relativity
                    .parameters
                    .internal
                    .radial_velocity
                    * (particle_normalized_position.y * element3s.z
                        - particle_normalized_position.z * element3s.y),
                3. * particle
                    .general_relativity
                    .parameters
                    .internal
                    .radial_velocity
                    * (particle_normalized_position.z * element3s.x
                        - particle_normalized_position.x * element3s.z),
                3. * particle
                    .general_relativity
                    .parameters
                    .internal
                    .radial_velocity
                    * (particle_normalized_position.x * element3s.y
                        - particle_normalized_position.y * element3s.x),
            );

            let factor_a = G / SPEED_OF_LIGHT_2;
            let mut total_general_relativity_spin_orbit_acceleration = Axes::from(
                factor_a * (element1.x - element2.x + element3.x),
                factor_a * (element1.y - element2.y + element3.y),
                factor_a * (element1.z - element2.z + element3.z),
            );

            particle
                .general_relativity
                .parameters
                .output
                .acceleration
                .add(&total_general_relativity_spin_orbit_acceleration);

            total_general_relativity_spin_orbit_acceleration
                .mul(particle.mass / host_particle.mass);
            sum_total_general_relativity_spin_orbit_acceleration
                .add(&total_general_relativity_spin_orbit_acceleration);

            //println!("{} {} {}", total_general_relativity_spin_orbit_acceleration.x, total_general_relativity_spin_orbit_acceleration.y, total_general_relativity_spin_orbit_acceleration.z);

            // Kidder 1995, equation 2.4a
            let mu = (host_particle.mass * particle.mass) / star_planet_mass;
            let mut newtonian_orbital_angular_momentum = particle
                .general_relativity
                .coordinates
                .position
                .cross(&particle.general_relativity.coordinates.velocity);
            newtonian_orbital_angular_momentum.mul(mu);

            let factor_mass = 2. + 3. / 2. * particle.mass / host_particle.mass;
            let mut element1 = newtonian_orbital_angular_momentum.cross(&star_angular_momentum);
            element1.mul(factor_mass);

            let element2 = particle_angular_momentum.cross(&star_angular_momentum);

            let scalar_product_particle_normalized_position_with_particle_angular_momentum =
                particle_normalized_position.dot(&particle_angular_momentum);

            let mut element3 = particle_normalized_position.cross(&star_angular_momentum);
            element3.mul(
                3. * scalar_product_particle_normalized_position_with_particle_angular_momentum,
            );

            host_particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .x += factor_a * (element1.x - element2.x + element3.x);
            host_particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .y += factor_a * (element1.y - element2.y + element3.y);
            host_particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .z += factor_a * (element1.z - element2.z + element3.z);
            //println!("{} {} {}", factor_a * (element1.x - element2.x + element3.x), factor_a * (element1.y - element2.y + element3.y), factor_a * (element1.z - element2.z + element3.z));

            // Kidder 1995, equation 2.4b
            let factor_mass = 2. + 3. / 2. * host_particle.mass / particle.mass;
            let mut element1 = newtonian_orbital_angular_momentum.cross(&particle_angular_momentum);
            element1.mul(factor_mass);

            let element2 = star_angular_momentum.cross(&particle_angular_momentum);

            let scalar_product_particle_normalized_position_with_star_angular_momentum =
                particle_normalized_position.dot(&star_angular_momentum);

            let mut element3 = particle_normalized_position.cross(&particle_angular_momentum);
            element3
                .mul(3. * scalar_product_particle_normalized_position_with_star_angular_momentum);

            particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .x = factor_a * (element1.x - element2.x + element3.x);
            particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .y = factor_a * (element1.y - element2.y + element3.y);
            particle
                .general_relativity
                .parameters
                .output
                .dangular_momentum_dt
                .z = factor_a * (element1.z - element2.z + element3.z);
        }
    }

    sum_total_general_relativity_spin_orbit_acceleration.negate();
    host_particle
        .general_relativity
        .parameters
        .output
        .acceleration
        .add(&sum_total_general_relativity_spin_orbit_acceleration);
}

////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
// [start] General Relativity based on REBOUNDx gr.c
pub fn calculate_anderson1975_general_relativity_acceleration(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    ignored_gravity_terms: IgnoreGravityTerms,
) {
    // Calculate Newtonian accelerations in the current setup and considering all particles

    let (host_newtonian_inertial_accelerations, newtonian_inertial_accelerations) =
        get_anderson1975_newhall1983_newtonian_inertial_accelerations(
            host_particle,
            particles,
            more_particles,
            ignored_gravity_terms,
        );

    // Transform to Jacobi coordinates
    let (
        jacobi_star_mass,
        _jacobi_star_position,
        _jacobi_star_velocity,
        _jacobi_star_acceleration,
        jacobi_particles_positions,
        jacobi_particles_velocities,
        mut jacobi_particles_accelerations,
    ) = anderson1975_general_relativity_inertial_to_jacobi_posvelacc(
        host_particle,
        particles,
        more_particles,
        host_newtonian_inertial_accelerations,
        newtonian_inertial_accelerations,
    );

    let n_particles = particles.len() + more_particles.len();
    let mu = host_particle.mass_g;
    for (
        jacobi_particle_acceleration,
        jacobi_particle_velocity,
        jacobi_particle_position,
        particle,
    ) in izip!(
        &mut jacobi_particles_accelerations[..n_particles],
        &jacobi_particles_velocities[..n_particles],
        &jacobi_particles_positions[..n_particles],
        particles.iter().chain(more_particles.iter())
    ) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let mut vi = *jacobi_particle_velocity;

            let mut vi2 = jacobi_particle_velocity.x.powi(2)
                + jacobi_particle_velocity.y.powi(2)
                + jacobi_particle_velocity.z.powi(2);
            let ri = sqrt!(
                jacobi_particle_position.x.powi(2)
                    + jacobi_particle_position.y.powi(2)
                    + jacobi_particle_position.z.powi(2)
            );
            let mut factor_a = (0.5 * vi2 + 3. * mu / ri) / SPEED_OF_LIGHT_2;
            let mut old_v;
            let mut dv;
            let max_iterations = 10;
            for q in 0..max_iterations {
                old_v = vi;
                vi = *jacobi_particle_velocity;
                vi.div(1. - factor_a);
                vi2 = vi.x * vi.x + vi.y * vi.y + vi.z * vi.z;
                factor_a = (0.5 * vi2 + 3. * mu / ri) / SPEED_OF_LIGHT_2;
                dv = vi;
                dv.sub(&old_v);
                if (dv.x * dv.x + dv.y * dv.y + dv.z * dv.z) / vi2 < DBL_EPSILON_2 {
                    break;
                } else if q == max_iterations {
                    println!(
                        "[WARNING {} UTC] {} iterations in general relativity failed to converge. This is typically because the perturbation is too strong for the current implementation.",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap(),
                        max_iterations
                    );
                }
            }

            let factor_b = (mu / ri - 1.5 * vi2) * mu / (ri * ri * ri) / SPEED_OF_LIGHT_2;
            let rdotrdot = jacobi_particle_position.dot(jacobi_particle_velocity);
            let vidot = Axes::from(
                jacobi_particle_acceleration.x + factor_b * jacobi_particle_position.x,
                jacobi_particle_acceleration.y + factor_b * jacobi_particle_position.y,
                jacobi_particle_acceleration.z + factor_b * jacobi_particle_position.z,
            );
            let vdotvdot = vi.dot(&vidot);
            let factor_d = (vdotvdot - 3. * mu / (ri * ri * ri) * rdotrdot) / SPEED_OF_LIGHT_2;
            jacobi_particle_acceleration.x =
                factor_b * (1. - factor_a) * jacobi_particle_position.x
                    - factor_a * jacobi_particle_acceleration.x
                    - factor_d * vi.x;
            jacobi_particle_acceleration.y =
                factor_b * (1. - factor_a) * jacobi_particle_position.y
                    - factor_a * jacobi_particle_acceleration.y
                    - factor_d * vi.y;
            jacobi_particle_acceleration.z =
                factor_b * (1. - factor_a) * jacobi_particle_position.z
                    - factor_a * jacobi_particle_acceleration.z
                    - factor_d * vi.z;
        }
    }

    let jacobi_star_acceleration = Axes::new();
    let (star_acceleration, particles_accelerations) =
        anderson1975_general_relativity_jacobi_to_inertial_acc(
            particles,
            more_particles,
            jacobi_star_mass,
            jacobi_star_acceleration,
            jacobi_particles_accelerations,
        );

    // This algorithm computes general_relativity.parameters.output.acceleration in the inertial frame,
    // which is the same coordinate system that is expressed all the rest of additional
    // effects
    for (particle, particle_acceleration) in particles
        .iter_mut()
        .chain(more_particles.iter_mut())
        .zip(particles_accelerations.iter())
    {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            particle.general_relativity.parameters.output.acceleration = *particle_acceleration;
        }
    }
    host_particle
        .general_relativity
        .parameters
        .output
        .acceleration = star_acceleration;
}

fn anderson1975_general_relativity_inertial_to_jacobi_posvelacc(
    host_particle: &Particle,
    particles: &[Particle],
    more_particles: &[Particle],
    host_newtonian_inertial_accelerations: Axes,
    newtonian_inertial_accelerations: [Axes; MAX_PARTICLES - 1],
) -> (
    f64,
    Axes,
    Axes,
    Axes,
    [Axes; MAX_PARTICLES - 1],
    [Axes; MAX_PARTICLES - 1],
    [Axes; MAX_PARTICLES - 1],
) {
    let mut jacobi_particles_positions = [Axes::new(); MAX_PARTICLES - 1];
    let mut jacobi_particles_velocities = [Axes::new(); MAX_PARTICLES - 1];
    let mut jacobi_particles_accelerations = [Axes::new(); MAX_PARTICLES - 1];

    let m0 = host_particle.mass;
    let mut eta = m0;
    let mut position = host_particle.inertial_position;
    position.mul(eta);
    let mut velocity = host_particle.inertial_velocity;
    velocity.mul(eta);
    let mut acceleration = host_newtonian_inertial_accelerations;
    acceleration.mul(eta);

    for (
        particle,
        particle_newtonian_inertial_accelerations,
        jacobi_particle_position,
        jacobi_particle_velocity,
        jacobi_particle_acceleration,
    ) in izip!(
        particles.iter().chain(more_particles.iter()), // zip will pick the lowest common number of elements
        &newtonian_inertial_accelerations,
        &mut jacobi_particles_positions,
        &mut jacobi_particles_velocities,
        &mut jacobi_particles_accelerations
    ) {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let ei = 1. / eta;
            eta += particle.mass;
            // multiply by ei first, then eta (pme = eta * ei);

            position.mul(ei);
            *jacobi_particle_position = particle.inertial_position;
            jacobi_particle_position.sub(&position);

            velocity.mul(ei);
            *jacobi_particle_velocity = particle.inertial_velocity;
            jacobi_particle_velocity.sub(&velocity);

            acceleration.mul(ei);
            *jacobi_particle_acceleration = *particle_newtonian_inertial_accelerations;
            jacobi_particle_acceleration.sub(&acceleration);

            position.mul(eta);
            let mut tmp = *jacobi_particle_position;
            tmp.mul(particle.mass);
            position.add(&tmp);

            velocity.mul(eta);
            let mut tmp = *jacobi_particle_velocity;
            tmp.mul(particle.mass);
            velocity.add(&tmp);

            acceleration.mul(eta);
            let mut tmp = *jacobi_particle_acceleration;
            tmp.mul(particle.mass);
            acceleration.add(&tmp);
        }
    }
    let jacobi_star_mass = eta;
    position.mul(1. / eta);
    velocity.mul(1. / eta);
    acceleration.mul(1. / eta);
    let jacobi_star_position = position;
    let jacobi_star_velocity = velocity;
    let jacobi_star_acceleration = acceleration;
    (
        jacobi_star_mass,
        jacobi_star_position,
        jacobi_star_velocity,
        jacobi_star_acceleration,
        jacobi_particles_positions,
        jacobi_particles_velocities,
        jacobi_particles_accelerations,
    )
}

fn anderson1975_general_relativity_jacobi_to_inertial_acc(
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    jacobi_star_mass: f64,
    jacobi_star_acceleration: Axes,
    jacobi_particles_accelerations: [Axes; MAX_PARTICLES - 1],
) -> (Axes, [Axes; MAX_PARTICLES - 1]) {
    let mut particles_accelerations = [Axes::new(); MAX_PARTICLES - 1];
    let n_particles = particles.len() + more_particles.len();

    let mut eta = jacobi_star_mass;
    let mut s = jacobi_star_acceleration;
    s.mul(eta);
    for ((particle, particle_acceleration), jacobi_particle_acceleration) in particles
        .iter()
        .chain(more_particles.iter())
        .rev()
        .zip(particles_accelerations[..n_particles].iter_mut().rev())
        .zip(jacobi_particles_accelerations[..n_particles].iter().rev())
    {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            let ei = 1. / eta;
            let mut tmp = *jacobi_particle_acceleration;
            tmp.mul(particle.mass);

            s.sub(&tmp);
            s.mul(ei);

            *particle_acceleration = *jacobi_particle_acceleration;
            particle_acceleration.add(&s);
            eta -= particle.mass;
            s.mul(eta);
        }
    }
    let mtot = eta;
    let mtot_i = 1. / mtot;
    let mut star_acceleration = s;
    star_acceleration.mul(mtot_i);
    (star_acceleration, particles_accelerations)
}
// [end] General Relativity based on REBOUNDx gr.c
//--------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////

fn get_anderson1975_newhall1983_newtonian_inertial_accelerations(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    ignored_gravity_terms: IgnoreGravityTerms,
) -> (Axes, [Axes; MAX_PARTICLES - 1]) {
    let mut newtonian_inertial_accelerations = [Axes::new(); MAX_PARTICLES - 1];
    let mut host_newtonian_inertial_accelerations = host_particle.inertial_acceleration;
    for (newtonian_acceleration, particle) in newtonian_inertial_accelerations
        .iter_mut()
        .zip(particles.iter_mut().chain(more_particles.iter_mut()))
    {
        *newtonian_acceleration = particle.inertial_acceleration;
    }
    // If some terms where ignored by the integrator, they should be added
    if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne
        || ignored_gravity_terms == IgnoreGravityTerms::WHFastTwo
    {
        let n_particles = if ignored_gravity_terms == IgnoreGravityTerms::WHFastOne {
            2 - 1 // Only host - next particle interaction, and the host is already out of the particles vector
        } else {
            particles.len() + more_particles.len()
        };
        for (newtonian_acceleration, particle) in newtonian_inertial_accelerations[..n_particles]
            .iter_mut()
            .zip(particles.iter_mut().chain(more_particles.iter_mut()))
        {
            let mut tmp = host_particle.inertial_position;
            tmp.sub(&particle.general_relativity.coordinates.position);
            let mut tmp2 = tmp;
            let r = tmp.norm();
            let prefac = G / r.powi(3);
            let prefac_mass_star = prefac * host_particle.mass;
            let prefac_mass_particle = prefac * particle.mass;

            tmp.mul(prefac_mass_particle);
            host_newtonian_inertial_accelerations.sub(&tmp);

            tmp2.mul(prefac_mass_star);
            newtonian_acceleration.add(&tmp2);
        }
    }
    (
        host_newtonian_inertial_accelerations,
        newtonian_inertial_accelerations,
    )
}

////////////////////////////////////////////////////////////////////////////
//--------------------------------------------------------------------------
// [start] General Relativity FULL based on REBOUNDx gr.c
pub fn calculate_newhall1983_general_relativity_acceleration(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    ignored_gravity_terms: IgnoreGravityTerms,
) {
    // host_particle is separated from particles for homogeneity with the rest of methods,
    // but this implementation of General Relativity does not uses a host, all
    // interactions between all the bodies are computed
    let mut host_a_const = [Axes::new()]; // array that stores the value of the constant term
    let mut host_a_new = [Axes::new()]; // stores the newly calculated term
    let mut host_rs = [[0.; MAX_PARTICLES]; 1];
    let mut host_drs = [[Axes::new(); MAX_PARTICLES]; 1];
    let mut a_const = [Axes::new(); MAX_PARTICLES]; // array that stores the value of the constant term
    let mut a_new = [Axes::new(); MAX_PARTICLES]; // stores the newly calculated term
    let mut rs = [[0.; MAX_PARTICLES]; MAX_PARTICLES];
    let mut drs = [[Axes::new(); MAX_PARTICLES]; MAX_PARTICLES];

    let (host_newtonian_inertial_accelerations, newtonian_inertial_accelerations) =
        get_anderson1975_newhall1983_newtonian_inertial_accelerations(
            host_particle,
            particles,
            more_particles,
            ignored_gravity_terms,
        );

    for (i, ((particle_i, drs_i), rs_i)) in iter::once(&*host_particle)
        .chain(particles.iter())
        .chain(more_particles.iter()) // zip will pick the lowest common number of elements
        .zip(host_drs.iter_mut().chain(drs.iter_mut()))
        .zip(host_rs.iter_mut().chain(rs.iter_mut()))
        .enumerate()
    {
        // compute distances
        for (j, (particle_j, drs_i_j)) in iter::once(&*host_particle)
            .chain(particles.iter())
            .chain(more_particles.iter()) // zip will pick the lowest common number of elements
            .zip(drs_i.iter_mut())
            .enumerate()
        {
            if j != i
                && (particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled
                    || particle_j.general_relativity.effect != GeneralRelativityEffect::Disabled)
            {
                *drs_i_j = particle_i.inertial_position;
                drs_i_j.sub(&particle_j.inertial_position);
                rs_i[j] = drs_i_j.norm();
                //println!("i j: {} {} {:e}", i, j, rs_i[j]);
            }
        }
    }

    for (i, (((particle_i, a_const_i), drs_i), rs_i)) in iter::once(&*host_particle)
        .chain(particles.iter())
        .chain(more_particles.iter()) // zip will pick the lowest common number of elements
        .zip(host_a_const.iter_mut().chain(a_const.iter_mut()))
        .zip(host_drs.iter_mut().chain(drs.iter_mut()))
        .zip(host_rs.iter().chain(rs.iter()))
        .enumerate()
    {
        // then compute the constant terms:
        let mut a_constx = 0.;
        let mut a_consty = 0.;
        let mut a_constz = 0.;
        // 1st constant part
        for (j, (particle_j, drs_i_j)) in iter::once(&*host_particle)
            .chain(particles.iter())
            .chain(more_particles.iter()) // zip will pick the lowest common number of elements
            .zip(drs_i.iter())
            .enumerate()
        {
            if j != i
                && (particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled
                    || particle_j.general_relativity.effect != GeneralRelativityEffect::Disabled)
            {
                let dxij = drs_i_j.x;
                let dyij = drs_i_j.y;
                let dzij = drs_i_j.z;
                let rij2 = rs_i[j].powi(2);
                let rij3 = rij2 * rs_i[j];

                let mut a1 = 0.;
                for (k, (particle_k, rs_i_k)) in iter::once(&*host_particle)
                    .chain(particles.iter())
                    .chain(more_particles.iter()) // zip will pick the lowest common number of elements
                    .zip(rs_i.iter())
                    .enumerate()
                {
                    if k != i {
                        a1 += (4. / (SPEED_OF_LIGHT_2)) * G * particle_k.mass / rs_i_k;
                    }
                }

                let mut a2 = 0.;
                for (l, (particle_l, rs_l)) in iter::once(&*host_particle)
                    .chain(particles.iter())
                    .chain(more_particles.iter()) // zip will pick the lowest common number of elements
                    .zip(host_rs.iter().chain(rs.iter()))
                    .enumerate()
                {
                    if l != j {
                        a2 += (1. / (SPEED_OF_LIGHT_2)) * G * particle_l.mass / rs_l[j];
                    }
                }

                let vi2 = particle_i.inertial_velocity.x.powi(2)
                    + particle_i.inertial_velocity.y.powi(2)
                    + particle_i.inertial_velocity.z.powi(2);
                let a3 = -vi2 / SPEED_OF_LIGHT_2;

                let vj2 = particle_j.inertial_velocity.x.powi(2)
                    + particle_j.inertial_velocity.y.powi(2)
                    + particle_j.inertial_velocity.z.powi(2);
                let a4 = -2. * vj2 / SPEED_OF_LIGHT_2;

                let a5 = (4. / SPEED_OF_LIGHT_2)
                    * (particle_i.inertial_velocity.x * particle_j.inertial_velocity.x
                        + particle_i.inertial_velocity.y * particle_j.inertial_velocity.y
                        + particle_i.inertial_velocity.z * particle_j.inertial_velocity.z);

                let a6_0 = dxij * particle_j.inertial_velocity.x
                    + dyij * particle_j.inertial_velocity.y
                    + dzij * particle_j.inertial_velocity.z;
                let a6 = (3. / (2. * SPEED_OF_LIGHT_2)) * a6_0.powi(2) / rij2;

                let factor1 = a1 + a2 + a3 + a4 + a5 + a6;
                //println!("factors {:e} {:e} {:e} {:e} {:e} {:e}", a1, a2, a3, a4, a5, a6);
                a_constx += G * particle_j.mass * dxij * factor1 / rij3;
                a_consty += G * particle_j.mass * dyij * factor1 / rij3;
                a_constz += G * particle_j.mass * dzij * factor1 / rij3;

                // 2nd constant part
                let dvxij = particle_i.inertial_velocity.x - particle_j.inertial_velocity.x;
                let dvyij = particle_i.inertial_velocity.y - particle_j.inertial_velocity.y;
                let dvzij = particle_i.inertial_velocity.z - particle_j.inertial_velocity.z;

                let factor2 = dxij
                    * (4. * particle_i.inertial_velocity.x - 3. * particle_j.inertial_velocity.x)
                    + dyij
                        * (4. * particle_i.inertial_velocity.y
                            - 3. * particle_j.inertial_velocity.y)
                    + dzij
                        * (4. * particle_i.inertial_velocity.z
                            - 3. * particle_j.inertial_velocity.z);

                a_constx += G * particle_j.mass * factor2 * dvxij / rij3 / SPEED_OF_LIGHT_2;
                a_consty += G * particle_j.mass * factor2 * dvyij / rij3 / SPEED_OF_LIGHT_2;
                a_constz += G * particle_j.mass * factor2 * dvzij / rij3 / SPEED_OF_LIGHT_2;
            }
        }
        a_const_i.x = a_constx;
        a_const_i.y = a_consty;
        a_const_i.z = a_constz;
        //println!("a_const_i {:?}", a_const_i);
    }

    let n_particles = particles.len() + more_particles.len();
    let dev_limit = 1.0e-30;
    let max_iterations = 10;
    // Now running the substitution again and again through the loop below
    for k in 0..max_iterations {
        let host_a_old = host_a_new;
        let a_old = a_new;
        // now add on the non-constant term
        for (i, ((((particle_i, a_new_i), drs_i), rs_i), a_const_i)) in iter::once(&*host_particle)
            .chain(particles.iter())
            .chain(more_particles.iter()) // zip will pick the lowest common number of elements
            .zip(host_a_new.iter_mut().chain(a_new.iter_mut()))
            .zip(host_drs.iter_mut().chain(drs.iter_mut()))
            .zip(host_rs.iter_mut().chain(rs.iter_mut()))
            .zip(host_a_const.iter().chain(a_const.iter()))
            .enumerate()
        {
            let mut non_constx = 0.;
            let mut non_consty = 0.;
            let mut non_constz = 0.;
            for (j, ((((particle_j, a_old_j), a_newton_j), drs_i_j), rs_i_j)) in
                iter::once(&*host_particle)
                    .chain(particles.iter())
                    .chain(more_particles.iter()) // zip will pick the lowest common number of elements
                    .zip(host_a_old.iter().chain(a_old.iter()))
                    .zip(
                        iter::once(&host_newtonian_inertial_accelerations)
                            .chain(newtonian_inertial_accelerations.iter()),
                    )
                    .zip(drs_i.iter())
                    .zip(rs_i.iter())
                    .enumerate()
            {
                if j != i
                    && (particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled
                        || particle_j.general_relativity.effect
                            != GeneralRelativityEffect::Disabled)
                {
                    let dxij = drs_i_j.x;
                    let dyij = drs_i_j.y;
                    let dzij = drs_i_j.z;
                    let rij = rs_i_j;
                    let rij2 = rij.powi(2);
                    let rij3 = rij2 * rij;
                    non_constx += (G * particle_j.mass * dxij / rij3)
                        * (dxij * (a_newton_j.x + a_old_j.x)
                            + dyij * (a_newton_j.y + a_old_j.y)
                            + dzij * (a_newton_j.z + a_old_j.z))
                        / (2. * SPEED_OF_LIGHT_2)
                        + (7. / (2. * SPEED_OF_LIGHT_2))
                            * G
                            * particle_j.mass
                            * (a_newton_j.x + a_old_j.x)
                            / rij;
                    non_consty += (G * particle_j.mass * dyij / rij3)
                        * (dxij * (a_newton_j.x + a_old_j.x)
                            + dyij * (a_newton_j.y + a_old_j.y)
                            + dzij * (a_newton_j.z + a_old_j.z))
                        / (2. * SPEED_OF_LIGHT_2)
                        + (7. / (2. * SPEED_OF_LIGHT_2))
                            * G
                            * particle_j.mass
                            * (a_newton_j.y + a_old_j.y)
                            / rij;
                    non_constz += (G * particle_j.mass * dzij / rij3)
                        * (dxij * (a_newton_j.x + a_old_j.x)
                            + dyij * (a_newton_j.y + a_old_j.y)
                            + dzij * (a_newton_j.z + a_old_j.z))
                        / (2. * SPEED_OF_LIGHT_2)
                        + (7. / (2. * SPEED_OF_LIGHT_2))
                            * G
                            * particle_j.mass
                            * (a_newton_j.z + a_old_j.z)
                            / rij;
                }
            }
            a_new_i.x = a_const_i.x + non_constx;
            a_new_i.y = a_const_i.y + non_consty;
            a_new_i.z = a_const_i.z + non_constz;
            //println!("non_constx {:?}", non_constx);
            //println!("non_consty {:?}", non_consty);
            //println!("non_constz {:?}", non_constz);
        }

        // break out loop if a_new is converging
        let mut maxdev = 0.;
        let mut dx = 0.;
        let mut dy = 0.;
        let mut dz = 0.;
        for ((particle_i, a_new_i), a_old_i) in iter::once(&*host_particle)
            .chain(particles.iter())
            .chain(more_particles.iter()) // zip will pick the lowest common number of elements
            .zip(host_a_new.iter_mut().chain(a_new[..n_particles].iter_mut()))
            .zip(host_a_old.iter().chain(a_old.iter()))
        {
            if particle_i.general_relativity.effect != GeneralRelativityEffect::Disabled {
                if abs!(a_new_i.x) < dev_limit {
                    dx = abs!(a_new_i.x - a_old_i.x) / a_new_i.x;
                }
                if abs!(a_new_i.y) < dev_limit {
                    dy = abs!(a_new_i.y - a_old_i.y) / a_new_i.y;
                }
                if abs!(a_new_i.z) < dev_limit {
                    dz = abs!(a_new_i.z - a_old_i.z) / a_new_i.z;
                }
                if dx > maxdev {
                    maxdev = dx;
                }
                if dy > maxdev {
                    maxdev = dy;
                }
                if dz > maxdev {
                    maxdev = dz;
                }
            }
        }

        if maxdev < dev_limit {
            break;
        } else if k == max_iterations {
            println!(
                "[WARNING {} UTC] {} iterations in general relativity failed to converge.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                max_iterations
            );
        }
    }

    //// update acceleration in particles
    // This algorithm computes general_relativity.parameters.output.acceleration in the inertial frame,
    // which is the same coordinate system that is expressed all the rest of additional
    // effects
    for (particle, a_new_particle) in particles
        .iter_mut()
        .chain(more_particles.iter_mut())
        .zip(a_new.iter())
    {
        if let GeneralRelativityEffect::OrbitingBody = particle.general_relativity.effect {
            particle.general_relativity.parameters.output.acceleration = *a_new_particle;
        }
    }
    host_particle
        .general_relativity
        .parameters
        .output
        .acceleration = host_a_new[0];
}
