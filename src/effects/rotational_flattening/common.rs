use super::super::super::Axes;
use super::super::super::Particle;
use super::super::tides::creep_coplanar::CreepCoplanarParameters;
use super::creep_coplanar;
use super::oblate_spheroid;
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleInternalParameters {
    pub distance: f64,
    // Oblate spheroid specific:
    pub scalar_product_of_vector_position_with_stellar_spin: f64,
    pub scalar_product_of_vector_position_with_planetary_spin: f64,
    pub radial_component_of_the_force_induced_by_rotation: f64,
    pub factor_for_the_force_induced_by_star_rotation: f64,
    pub factor_for_the_force_induced_by_planet_rotation: f64,
    pub orthogonal_component_of_the_force_induced_by_star_rotation: f64,
    pub orthogonal_component_of_the_force_induced_by_planet_rotation: f64,
    // Creep coplanar specific:
    pub shape: Axes,
    //
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleOutputParameters {
    pub acceleration: Axes,
    pub dangular_momentum_dt: Axes, // Force
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleParameters {
    pub internal: RotationalFlatteningParticleInternalParameters,
    pub output: RotationalFlatteningParticleOutputParameters,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlatteningParticleCoordinates {
    // Positions/velocities in a heliocentric frame
    // (i.e., the host is at rest with respect to the origin of the coordinate system)
    pub position: Axes,
    pub velocity: Axes,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum RotationalFlatteningModel {
    OblateSpheroid(oblate_spheroid::OblateSpheroidParameters),
    CreepCoplanar(CreepCoplanarParameters),
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum RotationalFlatteningEffect {
    CentralBody(RotationalFlatteningModel),
    OrbitingBody(RotationalFlatteningModel),
    Disabled,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct RotationalFlattening {
    pub effect: RotationalFlatteningEffect,
    pub parameters: RotationalFlatteningParticleParameters,
    pub coordinates: RotationalFlatteningParticleCoordinates,
}

impl RotationalFlattening {
    pub fn new(effect: RotationalFlatteningEffect) -> RotationalFlattening {
        RotationalFlattening {
            effect,
            parameters: RotationalFlatteningParticleParameters {
                internal: RotationalFlatteningParticleInternalParameters {
                    distance: 0.,
                    scalar_product_of_vector_position_with_stellar_spin: 0.,
                    scalar_product_of_vector_position_with_planetary_spin: 0.,
                    radial_component_of_the_force_induced_by_rotation: 0.,
                    factor_for_the_force_induced_by_star_rotation: 0.,
                    factor_for_the_force_induced_by_planet_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_star_rotation: 0.,
                    orthogonal_component_of_the_force_induced_by_planet_rotation: 0.,
                    shape: Axes::new(),
                },
                output: RotationalFlatteningParticleOutputParameters {
                    acceleration: Axes::new(),
                    dangular_momentum_dt: Axes::new(),
                },
            },
            coordinates: RotationalFlatteningParticleCoordinates {
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
    if let RotationalFlatteningEffect::CentralBody(_) = host_particle.rotational_flattening.effect {
        host_particle
            .rotational_flattening
            .parameters
            .internal
            .scalar_product_of_vector_position_with_stellar_spin = 0.;
        host_particle
            .rotational_flattening
            .parameters
            .internal
            .scalar_product_of_vector_position_with_planetary_spin = 0.;
        host_particle
            .rotational_flattening
            .parameters
            .output
            .acceleration
            .zero();
        host_particle
            .rotational_flattening
            .parameters
            .output
            .dangular_momentum_dt
            .zero();
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody(_) =
                particle.rotational_flattening.effect
            {
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_stellar_spin =
                    particle.rotational_flattening.coordinates.position.x * host_particle.spin.x
                        + particle.rotational_flattening.coordinates.position.y
                            * host_particle.spin.y
                        + particle.rotational_flattening.coordinates.position.z
                            * host_particle.spin.z;
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .scalar_product_of_vector_position_with_planetary_spin =
                    particle.rotational_flattening.coordinates.position.x * particle.spin.x
                        + particle.rotational_flattening.coordinates.position.y * particle.spin.y
                        + particle.rotational_flattening.coordinates.position.z * particle.spin.z;
                particle
                    .rotational_flattening
                    .parameters
                    .output
                    .acceleration
                    .zero();
                particle
                    .rotational_flattening
                    .parameters
                    .output
                    .dangular_momentum_dt
                    .zero();
            }
        }
    }
}

pub fn inertial_to_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let RotationalFlatteningEffect::CentralBody(_) = host_particle.rotational_flattening.effect {
        // Inertial to Heliocentric positions/velocities
        host_particle
            .rotational_flattening
            .coordinates
            .position
            .zero();
        host_particle
            .rotational_flattening
            .coordinates
            .velocity
            .zero();
        host_particle
            .rotational_flattening
            .parameters
            .internal
            .distance = 0.;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody(_) =
                particle.rotational_flattening.effect
            {
                particle.rotational_flattening.coordinates.position.x =
                    particle.inertial_position.x - host_particle.inertial_position.x;
                particle.rotational_flattening.coordinates.position.y =
                    particle.inertial_position.y - host_particle.inertial_position.y;
                particle.rotational_flattening.coordinates.position.z =
                    particle.inertial_position.z - host_particle.inertial_position.z;
                particle.rotational_flattening.coordinates.velocity.x =
                    particle.inertial_velocity.x - host_particle.inertial_velocity.x;
                particle.rotational_flattening.coordinates.velocity.y =
                    particle.inertial_velocity.y - host_particle.inertial_velocity.y;
                particle.rotational_flattening.coordinates.velocity.z =
                    particle.inertial_velocity.z - host_particle.inertial_velocity.z;
                particle.rotational_flattening.parameters.internal.distance = sqrt!(
                    particle
                        .rotational_flattening
                        .coordinates
                        .position
                        .x
                        .powi(2)
                        + particle
                            .rotational_flattening
                            .coordinates
                            .position
                            .y
                            .powi(2)
                        + particle
                            .rotational_flattening
                            .coordinates
                            .position
                            .z
                            .powi(2)
                );
            }
        }
    }
}

pub fn copy_heliocentric_coordinates(
    host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    if let RotationalFlatteningEffect::CentralBody(_) = host_particle.rotational_flattening.effect {
        host_particle
            .rotational_flattening
            .coordinates
            .position
            .copy_from(&host_particle.heliocentric_position);
        host_particle
            .rotational_flattening
            .coordinates
            .velocity
            .copy_from(&host_particle.heliocentric_velocity);
        host_particle
            .rotational_flattening
            .parameters
            .internal
            .distance = host_particle.heliocentric_distance;
        for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
            if let RotationalFlatteningEffect::OrbitingBody(_) =
                particle.rotational_flattening.effect
            {
                particle
                    .rotational_flattening
                    .coordinates
                    .position
                    .copy_from(&particle.heliocentric_position);
                particle
                    .rotational_flattening
                    .coordinates
                    .velocity
                    .copy_from(&particle.heliocentric_velocity);
                particle.rotational_flattening.parameters.internal.distance =
                    particle.heliocentric_distance;
            }
        }
    }
}

pub fn calculate_dangular_momentum_dt_induced_by_rotational_flattening(
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    let factor = -1.0;

    let central_body = false;
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::OrbitingBody(rotational_flattening_model) =
            particle.rotational_flattening.effect
        {
            let torque_induced_by_rotational_flattening = match rotational_flattening_model {
                RotationalFlatteningModel::OblateSpheroid(_) => {
                    oblate_spheroid::calculate_torque_induced_by_rotational_flattening(
                        rotational_flattening_host_particle,
                        particle,
                        central_body,
                    )
                }
                RotationalFlatteningModel::CreepCoplanar(_) => {
                    creep_coplanar::calculate_torque_induced_by_rotational_flattening(
                        rotational_flattening_host_particle,
                        particle,
                        central_body,
                    )
                }
            };
            // - Equation 25 from Bolmont et al. 2015
            // Integration of the spin (total torque rot):
            particle
                .rotational_flattening
                .parameters
                .output
                .dangular_momentum_dt
                .copy_from(&torque_induced_by_rotational_flattening);
            particle
                .rotational_flattening
                .parameters
                .output
                .dangular_momentum_dt
                .mul(factor);
        }
    }

    let central_body = true;
    let mut dangular_momentum_dt = Axes::new();
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::CentralBody(rotational_flattening_model) =
            rotational_flattening_host_particle
                .rotational_flattening
                .effect
        {
            let torque_induced_by_rotational_flattening = match rotational_flattening_model {
                RotationalFlatteningModel::OblateSpheroid(_) => {
                    oblate_spheroid::calculate_torque_induced_by_rotational_flattening(
                        rotational_flattening_host_particle,
                        particle,
                        central_body,
                    )
                }
                RotationalFlatteningModel::CreepCoplanar(_) => {
                    creep_coplanar::calculate_torque_induced_by_rotational_flattening(
                        rotational_flattening_host_particle,
                        particle,
                        central_body,
                    )
                }
            };
            // Integration of the spin (total torque rot):
            dangular_momentum_dt.x += factor * torque_induced_by_rotational_flattening.x;
            dangular_momentum_dt.y += factor * torque_induced_by_rotational_flattening.y;
            dangular_momentum_dt.z += factor * torque_induced_by_rotational_flattening.z;
        }
    }
    // - Equation 25 from Bolmont et al. 2015
    rotational_flattening_host_particle
        .rotational_flattening
        .parameters
        .output
        .dangular_momentum_dt
        .copy_from(&dangular_momentum_dt);
}

pub fn calculate_acceleration_induced_by_rotational_flattering(
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    let factor2 = 1. / rotational_flattening_host_particle.mass;

    // Calculation of the norm square of the spin for the star
    let mut sum_force_induced_by_rotation = Axes::new();

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::OrbitingBody(rotational_flattening_model) =
            particle.rotational_flattening.effect
        {
            let force_induced_by_rotation = match rotational_flattening_model {
                RotationalFlatteningModel::OblateSpheroid(_) => {
                    oblate_spheroid::calculate_acceleration_induced_by_rotational_flattering(
                        rotational_flattening_host_particle,
                        particle,
                    )
                }
                RotationalFlatteningModel::CreepCoplanar(_) => {
                    creep_coplanar::calculate_acceleration_induced_by_rotational_flattering(
                        rotational_flattening_host_particle,
                        particle,
                    )
                }
            };

            sum_force_induced_by_rotation.add(&force_induced_by_rotation);

            // - Equation 19 from Bolmont et al. 2015 (first term)
            let factor1 = 1. / particle.mass;
            particle
                .rotational_flattening
                .parameters
                .output
                .acceleration
                .copy_from(&force_induced_by_rotation);
            particle
                .rotational_flattening
                .parameters
                .output
                .acceleration
                .mul(factor1);
        }
    }

    // - Equation 19 from Bolmont et al. 2015 (second term)
    //for particle in particles.iter_mut() {
    //particle.rotational_flattening.parameters.output.acceleration.x += factor2 * sum_force_induced_by_rotation.x;
    //particle.rotational_flattening.parameters.output.acceleration.y += factor2 * sum_force_induced_by_rotation.y;
    //particle.rotational_flattening.parameters.output.acceleration.z += factor2 * sum_force_induced_by_rotation.z;
    //}
    // Instead of the previous code, keep star tidal acceleration separated:
    rotational_flattening_host_particle
        .rotational_flattening
        .parameters
        .output
        .acceleration
        .copy_from(&sum_force_induced_by_rotation);
    rotational_flattening_host_particle
        .rotational_flattening
        .parameters
        .output
        .acceleration
        .mul(-1.0 * factor2);
}
