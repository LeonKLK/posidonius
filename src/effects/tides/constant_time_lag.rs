use super::super::super::Axes;
use super::super::super::Particle;
use super::super::super::constants::{K2, MAX_PARTICLES, SMOOTHING_FACTOR_DYN_TIDE_COROTATION};
use super::super::super::tools;
use super::super::EvolutionType;
use super::common::TidalModel;
use super::common::TidesEffect;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct ConstantTimeLagParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64, // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the
                          // susceptibility of its shape to change in response to a tidal potential.
}

pub fn calculate_pair_dependent_scaled_dissipation_factors(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    calculate_planet_dependent_scaled_dissipation_factors(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
    calculate_host_dependent_scaled_dissipation_factors(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
}

fn calculate_planet_dependent_scaled_dissipation_factors(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    match tidal_host_particle.evolution {
        EvolutionType::BolmontMathis2016(_)
        | EvolutionType::GalletBolmont2017(_)
        | EvolutionType::LeconteChabrier2013(true) => {
            let host_norm_spin_vector = sqrt!(tidal_host_particle.norm_spin_vector_2);
            for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
                if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(_)) =
                    &particle.tides.effect
                {
                    let (
                        tidal_host_particle_dissipation_factor_scale,
                        tidal_host_particle_dissipation_factor,
                    ) = match &tidal_host_particle.tides.effect {
                        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(params)) => {
                            (params.dissipation_factor_scale, params.dissipation_factor)
                        }
                        _ => (0., 0.),
                    };
                    //
                    //// Excitation frequency needed by the model based on the
                    // instantaneous frequency (using positions, velocities and spins)
                    //let frequency = (particle.tides.coordinates.velocity.x - tidal_host_particle.spin.y*particle.tides.coordinates.position.z + tidal_host_particle.spin.z*particle.tides.coordinates.position.y).powi(2)
                    //+ (particle.tides.coordinates.velocity.y - tidal_host_particle.spin.z*particle.tides.coordinates.position.x + tidal_host_particle.spin.x*particle.tides.coordinates.position.z).powi(2)
                    //+ (particle.tides.coordinates.velocity.z - tidal_host_particle.spin.x*particle.tides.coordinates.position.y + tidal_host_particle.spin.y*particle.tides.coordinates.position.x).powi(2);
                    //let inverse_of_half_the_excitation_frequency = particle.tides.parameters.internal.distance / frequency;
                    // NOTE:  two_times_the_inverse_of_the_excitation_frequency: 2/w
                    //        inverse_of_half_the_excitation_frequency : 1/(w/2)
                    //
                    //// Excitation frequency needed by the model based on the
                    // mean frequency (using mean motion and spin).
                    //
                    // NOTE: The model is already here being used outside the
                    // validity domain, it seems not justified to use an
                    // instantaneous frequency.
                    let gm = tidal_host_particle.mass_g + particle.mass_g;
                    let (perihelion_distance, eccentricity) =
                        tools::calculate_perihelion_distance_and_eccentricity(
                            gm,
                            particle.tides.coordinates.position,
                            particle.tides.coordinates.velocity,
                        );
                    let mean_motion =
                        sqrt!(gm) * (perihelion_distance / (1.0 - eccentricity)).powf(-1.5);
                    let mut half_the_excitation_frequency =
                        abs!(host_norm_spin_vector - mean_motion);
                    // If the dynamical tide is excited, compute the planet-dependent stellar dissipation
                    // If the dynamical tide is not excited (i.e., equilibrium tide), do nothing since the default value corresponds to the planet scaled dissipation factor
                    if half_the_excitation_frequency < host_norm_spin_vector {
                        // The dynamical tide is excited
                        if half_the_excitation_frequency < SMOOTHING_FACTOR_DYN_TIDE_COROTATION {
                            half_the_excitation_frequency = SMOOTHING_FACTOR_DYN_TIDE_COROTATION;
                        }
                        let inverse_of_half_the_excitation_frequency =
                            1. / half_the_excitation_frequency;

                        // Eq. 4 and Eq.10 of Bolmont & Mathis have a typo, see page 5 of Gallet &
                        // Bolmont 2017. The lag angle has a 1/k2, but here we should *k2, that is why
                        // k2 does not appear here
                        // We add here the dissipation from the equilibrium tide to the dynamical tide one
                        let particle_dependent_scaled_dissipation_factor =
                            tidal_host_particle_dissipation_factor_scale
                                * (2.0 * K2 / (3.0 * tidal_host_particle.radius.powi(5))
                                    * tidal_host_particle.tides.parameters.internal.lag_angle
                                    * inverse_of_half_the_excitation_frequency
                                    + tidal_host_particle_dissipation_factor);

                        set_pair_dependent_scaled_dissipation_factor(
                            pair_dependent_scaled_dissipation_factor,
                            tidal_host_particle.id,
                            particle.id,
                            particle_dependent_scaled_dissipation_factor,
                        );
                    } else {
                        // Equilibrium regime (not dynamical tide)
                        remove_pair_dependent_scaled_dissipation_factor(
                            pair_dependent_scaled_dissipation_factor,
                            tidal_host_particle.id,
                            particle.id,
                        );
                    }
                }
            }
            //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
        }
        _ => {}
    }
}

fn calculate_host_dependent_scaled_dissipation_factors(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        match particle.evolution {
            EvolutionType::BolmontMathis2016(_)
            | EvolutionType::GalletBolmont2017(_)
            | EvolutionType::LeconteChabrier2013(true) => {
                if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(params)) =
                    &particle.tides.effect
                {
                    let particle_norm_spin_vector = sqrt!(particle.norm_spin_vector_2);
                    //// Excitation frequency needed by the model based on the
                    // mean frequency (using mean motion and spin).
                    //
                    // NOTE: The model is already here being used outside the
                    // validity domain, it seems not justified to use an
                    // instantaneous frequency.
                    let gm = tidal_host_particle.mass_g + particle.mass_g;
                    let (perihelion_distance, eccentricity) =
                        tools::calculate_perihelion_distance_and_eccentricity(
                            gm,
                            particle.tides.coordinates.position,
                            particle.tides.coordinates.velocity,
                        );
                    let mean_motion =
                        sqrt!(gm) * (perihelion_distance / (1.0 - eccentricity)).powf(-1.5);
                    let mut half_the_excitation_frequency =
                        abs!(particle_norm_spin_vector - mean_motion);
                    // If the dynamical tide is excited, compute the planet-dependent stellar dissipation
                    // If the dynamical tide is not excited (i.e., equilibrium tide), do nothing since the default value corresponds to the planet scaled dissipation factor
                    if half_the_excitation_frequency < particle_norm_spin_vector {
                        // The dynamical tide is excited
                        if half_the_excitation_frequency < SMOOTHING_FACTOR_DYN_TIDE_COROTATION {
                            half_the_excitation_frequency = SMOOTHING_FACTOR_DYN_TIDE_COROTATION;
                        }
                        let inverse_of_half_the_excitation_frequency =
                            1. / half_the_excitation_frequency;

                        // Eq. 4 and Eq.10 of Bolmont & Mathis have a typo, see page 5 of Gallet &
                        // Bolmont 2017. The lag angle has a 1/k2, but here we should *k2, that is why
                        // k2 does not appear here
                        // We add here the dissipation from the equilibrium tide to the dynamical tide one
                        let host_dependent_scaled_dissipation_factor = params
                            .dissipation_factor_scale
                            * (2.0 * K2 / (3.0 * particle.radius.powi(5))
                                * particle.tides.parameters.internal.lag_angle
                                * inverse_of_half_the_excitation_frequency
                                + params.dissipation_factor);

                        set_pair_dependent_scaled_dissipation_factor(
                            pair_dependent_scaled_dissipation_factor,
                            particle.id,
                            tidal_host_particle.id,
                            host_dependent_scaled_dissipation_factor,
                        );
                    } else {
                        // Equilibrium regime (not dynamical tide)
                        remove_pair_dependent_scaled_dissipation_factor(
                            pair_dependent_scaled_dissipation_factor,
                            tidal_host_particle.id,
                            particle.id,
                        );
                    }
                }
            }
            _ => {}
        }
    }
    //panic!("Please, contact Posidonius authors before using BolmontMathis2016/GalletBolmont2017/LeconteChabrier2013(true) evolutionary models. They may not be ready yet for scientific explotation.")
}

pub fn set_pair_dependent_scaled_dissipation_factor(
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
    id: usize,
    depends_on_id: usize,
    scaled_dissipation_factor: f64,
) {
    let key = id * MAX_PARTICLES + depends_on_id;
    pair_dependent_scaled_dissipation_factor.insert(key, scaled_dissipation_factor);
}

pub fn remove_pair_dependent_scaled_dissipation_factor(
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
    id: usize,
    depends_on_id: usize,
) {
    let key = id * MAX_PARTICLES + depends_on_id;
    pair_dependent_scaled_dissipation_factor.remove(&key);
}

pub fn get_pair_dependent_scaled_dissipation_factor_or_else(
    pair_dependent_scaled_dissipation_factor: &HashMap<usize, f64>,
    id: usize,
    depends_on_id: usize,
    evolution: EvolutionType,
    scaled_dissipation_factor: f64,
) -> f64 {
    let key = id * MAX_PARTICLES + depends_on_id;
    match evolution {
        EvolutionType::BolmontMathis2016(_)
        | EvolutionType::GalletBolmont2017(_)
        | EvolutionType::LeconteChabrier2013(true) => {
            match pair_dependent_scaled_dissipation_factor.get(&key) {
                Some(&value) => value,
                _ => scaled_dissipation_factor, // This happens if it goes out of the dynamical regime into the equilibrium regime
            }
        }
        _ => scaled_dissipation_factor,
    }
}

//////////////////////////////////////////////////////////////////////////////
/// TIDES
pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    let mut reference_spin = Axes::new();
    let orthogonal_component_of_the_tidal_force;
    let reference_rscalspin;

    if central_body {
        reference_spin.copy_from(&tidal_host_particle.spin);
        reference_rscalspin = particle
            .tides
            .parameters
            .internal
            .scalar_product_of_vector_position_with_stellar_spin;
        orthogonal_component_of_the_tidal_force = particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide;
    } else {
        reference_spin.copy_from(&particle.spin);
        reference_rscalspin = particle
            .tides
            .parameters
            .internal
            .scalar_product_of_vector_position_with_planetary_spin;
        orthogonal_component_of_the_tidal_force = particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide;
    }

    // distance to star
    let distance = particle.tides.parameters.internal.distance;

    //// Torque calculation (star)
    // - Equation 8-9 from Bolmont et al. 2015
    let torque_due_to_tides_x = orthogonal_component_of_the_tidal_force
        * (distance * reference_spin.x
            - reference_rscalspin * particle.tides.coordinates.position.x / distance
            - 1.0 / distance
                * (particle.tides.coordinates.position.y * particle.tides.coordinates.velocity.z
                    - particle.tides.coordinates.position.z
                        * particle.tides.coordinates.velocity.y));

    let torque_due_to_tides_y = orthogonal_component_of_the_tidal_force
        * (distance * reference_spin.y
            - reference_rscalspin * particle.tides.coordinates.position.y / distance
            - 1.0 / distance
                * (particle.tides.coordinates.position.z * particle.tides.coordinates.velocity.x
                    - particle.tides.coordinates.position.x
                        * particle.tides.coordinates.velocity.z));

    let torque_due_to_tides_z = orthogonal_component_of_the_tidal_force
        * (distance * reference_spin.z
            - reference_rscalspin * particle.tides.coordinates.position.z / distance
            - 1.0 / distance
                * (particle.tides.coordinates.position.x * particle.tides.coordinates.velocity.y
                    - particle.tides.coordinates.position.y
                        * particle.tides.coordinates.velocity.x));

    Axes::from(
        torque_due_to_tides_x,
        torque_due_to_tides_y,
        torque_due_to_tides_z,
    )
}

pub fn calculate_orthogonal_component_of_the_tidal_force(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    let central_body = true;
    calculate_orthogonal_component_of_the_tidal_force_for(
        central_body,
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
    calculate_orthogonal_component_of_the_tidal_force_for(
        !central_body,
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
}

fn calculate_orthogonal_component_of_the_tidal_force_for(
    central_body: bool,
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(_)) = &particle.tides.effect {
            //// Only calculate tides if planet is not in disk
            //if particle.disk_interaction_time == 0.0 {

            // (distance to star)^7
            let distance_7 = particle.tides.parameters.internal.distance.powi(7);

            //// Tidal force calculation (star) :: Only orthogonal component is needed
            if central_body {
                // - Third line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5)
                //   because it uses sigma (i.e., scaled_dissipation_factor)
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                let host_scaled_dissipation_factor =
                    get_pair_dependent_scaled_dissipation_factor_or_else(
                        pair_dependent_scaled_dissipation_factor,
                        tidal_host_particle.id,
                        particle.id,
                        tidal_host_particle.evolution,
                        tidal_host_particle
                            .tides
                            .parameters
                            .internal
                            .scaled_dissipation_factor,
                    );
                particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 4.5
                    * (particle.mass.powi(2))
                    * (tidal_host_particle.radius.powi(10))
                    * host_scaled_dissipation_factor
                    / distance_7;
            } else {
                // - Second line of Equation 5 from Bolmont et al. 2015
                //   This expression has R**10 (instead of R**5 in Eq. 5)
                //   because it uses sigma (i.e., scaled_dissipation_factor)
                //   and not k2$\Delta$t (between k2$\Delta$t and sigma
                //   there is a R**5 factor as shown in Equation 28)
                //   - k2 is love number
                let particle_scaled_dissipation_factor =
                    get_pair_dependent_scaled_dissipation_factor_or_else(
                        pair_dependent_scaled_dissipation_factor,
                        particle.id,
                        tidal_host_particle.id,
                        particle.evolution,
                        particle.tides.parameters.internal.scaled_dissipation_factor,
                    );
                particle
                    .tides
                    .parameters
                    .internal
                    .orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 4.5
                    * (tidal_host_particle.mass.powi(2))
                    * (particle.radius.powi(10))
                    * particle_scaled_dissipation_factor
                    / distance_7;

                // SBC
                //println!("> {:e} {:e} {:e} {:e}", tidal_host_particle.mass_g, particle.radius.powi(10), particle.scaled_dissipation_factor, distance_7);
                //println!("> {:e} {:e} {:e}", particle.tides.coordinates.position.x, particle.tides.coordinates.position.y, particle.tides.coordinates.position.z);
            }
            //} else {
            //if central_body {
            //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide = 0.0
            //} else{
            //particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide = 0.0
            //}
            //}
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    let host_mass_2 = tidal_host_particle.mass * tidal_host_particle.mass;

    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(params)) =
            &particle.tides.effect
        {
            let tidal_host_particle_love_number = match &tidal_host_particle.tides.effect {
                TidesEffect::CentralBody(TidalModel::ConstantTimeLag(params)) => params.love_number,
                _ => 0.,
            };
            let particle_mass_2 = particle.mass * particle.mass;
            // Conservative part of the radial tidal force
            let radial_component_of_the_tidal_force_conservative_part = -3.0 * K2
                / particle.tides.parameters.internal.distance.powi(7)
                * (particle_mass_2
                    * tidal_host_particle.radius.powi(5)
                    * tidal_host_particle_love_number
                    + host_mass_2 * particle.radius.powi(5) * params.love_number);

            // Dissipative part of the radial tidal force:
            let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity
                / particle.tides.parameters.internal.distance.powi(8);
            let host_scaled_dissipation_factor =
                get_pair_dependent_scaled_dissipation_factor_or_else(
                    pair_dependent_scaled_dissipation_factor,
                    tidal_host_particle.id,
                    particle.id,
                    tidal_host_particle.evolution,
                    tidal_host_particle
                        .tides
                        .parameters
                        .internal
                        .scaled_dissipation_factor,
                );
            let particle_scaled_dissipation_factor =
                get_pair_dependent_scaled_dissipation_factor_or_else(
                    pair_dependent_scaled_dissipation_factor,
                    particle.id,
                    tidal_host_particle.id,
                    particle.evolution,
                    particle.tides.parameters.internal.scaled_dissipation_factor,
                );
            let term1 = particle_mass_2
                * tidal_host_particle.radius.powi(10)
                * host_scaled_dissipation_factor;
            let term2 = host_mass_2 * particle.radius.powi(10) * particle_scaled_dissipation_factor;
            // If we consider the star as a point mass (used for denergy_dt calculation):
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass =
                factor1 * term2;
            let radial_component_of_the_tidal_force_dissipative_part = particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass
                + factor1 * term1;

            // Sum of the dissipative and conservative part of the radial force
            // - First line Equation 5 from Bolmont et al. 2015
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force =
                radial_component_of_the_tidal_force_conservative_part
                    + radial_component_of_the_tidal_force_dissipative_part;
        }
    }
}

pub fn calculate_tidal_force(tidal_host_particle: &Particle, particle: &Particle) -> Axes {
    // - Equation 6 from Bolmont et al. 2015
    let factor3 = particle
        .tides
        .parameters
        .internal
        .radial_component_of_the_tidal_force
        + (particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
            + particle
                .tides
                .parameters
                .internal
                .orthogonal_component_of_the_tidal_force_due_to_planetary_tide)
            * particle.tides.parameters.internal.radial_velocity
            / particle.tides.parameters.internal.distance;
    let total_tidal_force_x = factor3 * particle.tides.coordinates.position.x
        / particle.tides.parameters.internal.distance
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
            / particle.tides.parameters.internal.distance
            * (tidal_host_particle.spin.y * particle.tides.coordinates.position.z
                - tidal_host_particle.spin.z * particle.tides.coordinates.position.y
                - particle.tides.coordinates.velocity.x)
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
            / particle.tides.parameters.internal.distance
            * (particle.spin.y * particle.tides.coordinates.position.z
                - particle.spin.z * particle.tides.coordinates.position.y
                - particle.tides.coordinates.velocity.x);
    let total_tidal_force_y = factor3 * particle.tides.coordinates.position.y
        / particle.tides.parameters.internal.distance
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
            / particle.tides.parameters.internal.distance
            * (tidal_host_particle.spin.z * particle.tides.coordinates.position.x
                - tidal_host_particle.spin.x * particle.tides.coordinates.position.z
                - particle.tides.coordinates.velocity.y)
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
            / particle.tides.parameters.internal.distance
            * (particle.spin.z * particle.tides.coordinates.position.x
                - particle.spin.x * particle.tides.coordinates.position.z
                - particle.tides.coordinates.velocity.y);
    let total_tidal_force_z = factor3 * particle.tides.coordinates.position.z
        / particle.tides.parameters.internal.distance
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide
            / particle.tides.parameters.internal.distance
            * (tidal_host_particle.spin.x * particle.tides.coordinates.position.y
                - tidal_host_particle.spin.y * particle.tides.coordinates.position.x
                - particle.tides.coordinates.velocity.z)
        + particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide
            / particle.tides.parameters.internal.distance
            * (particle.spin.x * particle.tides.coordinates.position.y
                - particle.spin.y * particle.tides.coordinates.position.x
                - particle.tides.coordinates.velocity.z);
    //println!("factor3 {:e} {:e} {:e}", particle.tides.parameters.internal.radial_component_of_the_tidal_force, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_stellar_tide, particle.tides.parameters.internal.orthogonal_component_of_the_tidal_force_due_to_planetary_tide);
    //println!("d {:e} vrad {:e}", particle.tides.parameters.internal., particle.tides.parameters.internal.radial_velocity);
    //println!("total {:e} {:e} {:e}", total_tidal_force_x, total_tidal_force_y, total_tidal_force_z);
    Axes::from(
        total_tidal_force_x,
        total_tidal_force_y,
        total_tidal_force_z,
    )
}
