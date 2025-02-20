use super::super::super::Axes;
use super::super::super::Particle;
use super::common::RotationalFlatteningEffect;
use super::common::RotationalFlatteningModel;
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct OblateSpheroidParameters {
    pub love_number: f64, // love number for a completely fluid planet (used for rotational flattening effects)
}

pub fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening(
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    let central_body = true;
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(
        central_body,
        rotational_flattening_host_particle,
        particles,
        more_particles,
    );
    calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(
        !central_body,
        rotational_flattening_host_particle,
        particles,
        more_particles,
    );
}

fn calculate_orthogonal_component_of_the_force_induced_by_rotational_flattening_for(
    central_body: bool,
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    // Calculation of the norm square of the spin for the planet
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::OrbitingBody(
            RotationalFlatteningModel::OblateSpheroid(params),
        ) = particle.rotational_flattening.effect
        {
            let rotational_flattening_host_particle_love_number =
                match rotational_flattening_host_particle
                    .rotational_flattening
                    .effect
                {
                    RotationalFlatteningEffect::CentralBody(
                        RotationalFlatteningModel::OblateSpheroid(params),
                    ) => params.love_number,
                    _ => 0.,
                };
            if central_body {
                // - Star Equation 16 from Bolmont et al. 2015
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .factor_for_the_force_induced_by_star_rotation = particle.mass
                    * rotational_flattening_host_particle_love_number
                    * rotational_flattening_host_particle.norm_spin_vector_2
                    * rotational_flattening_host_particle.radius.powi(5)
                    / 6.; // Msun.AU^5.day-2
                // - Second part of Equation 15 from Bolmont et al. 2015
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .orthogonal_component_of_the_force_induced_by_star_rotation = -6.
                    * particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .factor_for_the_force_induced_by_star_rotation
                    * particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .scalar_product_of_vector_position_with_stellar_spin
                    / (rotational_flattening_host_particle.norm_spin_vector_2
                        * particle
                            .rotational_flattening
                            .parameters
                            .internal
                            .distance
                            .powi(5));
            } else {
                // - Planet Equation 16 from Bolmont et al. 2015
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .factor_for_the_force_induced_by_planet_rotation =
                    rotational_flattening_host_particle.mass
                        * params.love_number
                        * particle.norm_spin_vector_2
                        * particle.radius.powi(5)
                        / 6.; // Msun.AU^5.day-2
                // - Second part of Equation 15 from Bolmont et al. 2015
                particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .orthogonal_component_of_the_force_induced_by_planet_rotation = -6.
                    * particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .factor_for_the_force_induced_by_planet_rotation
                    * particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .scalar_product_of_vector_position_with_planetary_spin
                    / (particle.norm_spin_vector_2
                        * particle
                            .rotational_flattening
                            .parameters
                            .internal
                            .distance
                            .powi(5));
            }
        }
    }
}

pub fn calculate_radial_component_of_the_force_induced_by_rotational_flattening(
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let RotationalFlatteningEffect::OrbitingBody(_) = particle.rotational_flattening.effect {
            // - First part of Equation 15 from Bolmont et al. 2015
            particle
                .rotational_flattening
                .parameters
                .internal
                .radial_component_of_the_force_induced_by_rotation = -3.
                / particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .distance
                    .powi(5)
                * (particle
                    .rotational_flattening
                    .parameters
                    .internal
                    .factor_for_the_force_induced_by_planet_rotation
                    + particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .factor_for_the_force_induced_by_star_rotation)
                + 15.
                    / particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .distance
                        .powi(7)
                    * (particle
                        .rotational_flattening
                        .parameters
                        .internal
                        .factor_for_the_force_induced_by_star_rotation
                        * particle
                            .rotational_flattening
                            .parameters
                            .internal
                            .scalar_product_of_vector_position_with_stellar_spin
                        * particle
                            .rotational_flattening
                            .parameters
                            .internal
                            .scalar_product_of_vector_position_with_stellar_spin
                        / rotational_flattening_host_particle.norm_spin_vector_2
                        + particle
                            .rotational_flattening
                            .parameters
                            .internal
                            .factor_for_the_force_induced_by_planet_rotation
                            * particle
                                .rotational_flattening
                                .parameters
                                .internal
                                .scalar_product_of_vector_position_with_planetary_spin
                            * particle
                                .rotational_flattening
                                .parameters
                                .internal
                                .scalar_product_of_vector_position_with_planetary_spin
                            / particle.norm_spin_vector_2); // Msun.AU.day-1
        }
    }
}

pub fn calculate_torque_induced_by_rotational_flattening(
    rotational_flattening_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    let mut reference_spin = rotational_flattening_host_particle.spin;
    let orthogonal_component_of_the_force_induced_by_rotation;

    if central_body {
        orthogonal_component_of_the_force_induced_by_rotation = particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_star_rotation;
    } else {
        reference_spin = particle.spin;
        orthogonal_component_of_the_force_induced_by_rotation = particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_planet_rotation;
    }

    //// Torque calculation due to rotational flattening
    // - Equation 17-18 from Bolmont et al. 2015
    let torque_induced_by_rotational_flattening_x =
        orthogonal_component_of_the_force_induced_by_rotation
            * (particle.rotational_flattening.coordinates.position.y * reference_spin.z
                - particle.rotational_flattening.coordinates.position.z * reference_spin.y);
    let torque_induced_by_rotational_flattening_y =
        orthogonal_component_of_the_force_induced_by_rotation
            * (particle.rotational_flattening.coordinates.position.z * reference_spin.x
                - particle.rotational_flattening.coordinates.position.x * reference_spin.z);
    let torque_induced_by_rotational_flattening_z =
        orthogonal_component_of_the_force_induced_by_rotation
            * (particle.rotational_flattening.coordinates.position.x * reference_spin.y
                - particle.rotational_flattening.coordinates.position.y * reference_spin.x);

    Axes::from(
        torque_induced_by_rotational_flattening_x,
        torque_induced_by_rotational_flattening_y,
        torque_induced_by_rotational_flattening_z,
    )
}

pub fn calculate_acceleration_induced_by_rotational_flattering(
    rotational_flattening_host_particle: &Particle,
    particle: &Particle,
) -> Axes {
    // - Equation 15 from Bolmont et al. 2015
    let total_force_induced_by_rotation_x = particle
        .rotational_flattening
        .parameters
        .internal
        .radial_component_of_the_force_induced_by_rotation
        * particle.rotational_flattening.coordinates.position.x
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_planet_rotation
            * particle.spin.x
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_star_rotation
            * rotational_flattening_host_particle.spin.x;
    let total_force_induced_by_rotation_y = particle
        .rotational_flattening
        .parameters
        .internal
        .radial_component_of_the_force_induced_by_rotation
        * particle.rotational_flattening.coordinates.position.y
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_planet_rotation
            * particle.spin.y
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_star_rotation
            * rotational_flattening_host_particle.spin.y;
    let total_force_induced_by_rotation_z = particle
        .rotational_flattening
        .parameters
        .internal
        .radial_component_of_the_force_induced_by_rotation
        * particle.rotational_flattening.coordinates.position.z
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_planet_rotation
            * particle.spin.z
        + particle
            .rotational_flattening
            .parameters
            .internal
            .orthogonal_component_of_the_force_induced_by_star_rotation
            * rotational_flattening_host_particle.spin.z;
    //println!("Rot force = {:e} {:e} {:e}", total_force_induced_by_rotation_x, total_force_induced_by_rotation_y, total_force_induced_by_rotation_z);

    Axes::from(
        total_force_induced_by_rotation_x,
        total_force_induced_by_rotation_y,
        total_force_induced_by_rotation_z,
    )
}
