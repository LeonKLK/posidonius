use super::{TidalModel, TidesEffect};
use crate::EvolutionType;
use crate::Particle;
use crate::constants::{K2, MAX_PARTICLES, SMOOTHING_FACTOR_DYN_TIDE_COROTATION};
use crate::particles::Axes;
use crate::tools;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Selects which dissipation contributions the CTL model applies for a body.
///
/// The equilibrium tide dissipates through the constant `dissipation_factor`
/// (sigma). The frequency-averaged dynamical tide (Bolmont & Mathis 2016)
/// dissipates through the evolving `lag_angle` (derived from the inverse
/// tidal quality factor of the stellar evolution tables) and is only active
/// inside the excitation regime `|spin - n| < spin`; it therefore requires an
/// evolution model providing 1/Q (BolmontMathis2016, GalletBolmont2017 or
/// LeconteChabrier2013 with dissipation of dynamical tides).
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
pub enum TideComposition {
    /// Only the equilibrium tide, regardless of the evolution model.
    Equilibrium,
    /// Only the frequency-averaged dynamical tide; zero dissipation outside
    /// the excitation regime or when the evolution model provides no 1/Q.
    Dynamical,
    /// Dynamical tide added on top of the equilibrium tide when excited,
    /// equilibrium tide otherwise (historical behavior, default).
    #[default]
    Both,
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct ConstantTimeLagParameters {
    pub dissipation_factor: f64,
    pub dissipation_factor_scale: f64, // to scale the dissipation factor (multiply)
    pub love_number: f64, // Love number of degree 2 (i.e., k2). Dimensionless parameters that measure the rigidity of a planetary body and the
    // susceptibility of its shape to change in response to a tidal potential.
    #[serde(default)]
    pub tide_composition: TideComposition,
}

/// True for the evolution models that provide the frequency-averaged inverse
/// tidal quality factor needed by the dynamical tide (through `lag_angle`).
fn evolution_provides_inverse_tidal_q(evolution: EvolutionType) -> bool {
    matches!(
        evolution,
        EvolutionType::BolmontMathis2016(_)
            | EvolutionType::GalletBolmont2017(_)
            | EvolutionType::LeconteChabrier2013(true)
    )
}

/// Sigma contributed by the equilibrium tide alone.
fn equilibrium_scaled_dissipation_factor(params: &ConstantTimeLagParameters) -> f64 {
    params.dissipation_factor_scale * params.dissipation_factor
}

/// Sigma contributed by the frequency-averaged dynamical tide alone.
///
/// Eq. 4 and Eq. 10 of Bolmont & Mathis 2016 have a typo, see page 5 of
/// Gallet & Bolmont 2017. The lag angle has a 1/k2, but here we should *k2,
/// that is why k2 does not appear here.
fn dynamical_scaled_dissipation_factor(
    params: &ConstantTimeLagParameters,
    radius: f64,
    lag_angle: f64,
    inverse_of_half_the_excitation_frequency: f64,
) -> f64 {
    params.dissipation_factor_scale
        * (2.0 * K2 / (3.0 * radius.powi(5)))
        * lag_angle
        * inverse_of_half_the_excitation_frequency
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
    // Dynamical-tide sigma overrides only exist for a CTL central body whose
    // tide composition includes the dynamical tide and whose evolution model
    // provides the frequency-averaged 1/Q (through lag_angle).
    let host_params = match &tidal_host_particle.tides.effect {
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(params)) => *params,
        _ => return,
    };
    if host_params.tide_composition == TideComposition::Equilibrium
        || !evolution_provides_inverse_tidal_q(tidal_host_particle.evolution)
    {
        return;
    }
    let host_norm_spin_vector = sqrt!(tidal_host_particle.norm_spin_vector_2);
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(_)) = &particle.tides.effect {
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
            let mean_motion = sqrt!(gm) * (perihelion_distance / (1.0 - eccentricity)).powf(-1.5);
            let mut half_the_excitation_frequency = abs!(host_norm_spin_vector - mean_motion);
            // If the dynamical tide is excited, compute the planet-dependent stellar dissipation
            // If the dynamical tide is not excited (i.e., equilibrium tide), do nothing since the default value corresponds to the planet scaled dissipation factor
            if half_the_excitation_frequency < host_norm_spin_vector {
                // The dynamical tide is excited
                if half_the_excitation_frequency < SMOOTHING_FACTOR_DYN_TIDE_COROTATION {
                    half_the_excitation_frequency = SMOOTHING_FACTOR_DYN_TIDE_COROTATION;
                }
                let inverse_of_half_the_excitation_frequency = 1. / half_the_excitation_frequency;

                let dynamical = dynamical_scaled_dissipation_factor(
                    &host_params,
                    tidal_host_particle.radius,
                    tidal_host_particle.tides.parameters.internal.lag_angle,
                    inverse_of_half_the_excitation_frequency,
                );
                let particle_dependent_scaled_dissipation_factor =
                    match host_params.tide_composition {
                        TideComposition::Dynamical => dynamical,
                        // Historical behavior: dissipation from the equilibrium
                        // tide is added to the dynamical tide one
                        _ => dynamical + equilibrium_scaled_dissipation_factor(&host_params),
                    };

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

fn calculate_host_dependent_scaled_dissipation_factors(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // Dynamical-tide sigma overrides only exist for a CTL orbiting body
        // whose tide composition includes the dynamical tide and whose
        // evolution model provides the frequency-averaged 1/Q.
        let params = match &particle.tides.effect {
            TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(params)) => *params,
            _ => continue,
        };
        if params.tide_composition == TideComposition::Equilibrium
            || !evolution_provides_inverse_tidal_q(particle.evolution)
        {
            continue;
        }
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
        let mean_motion = sqrt!(gm) * (perihelion_distance / (1.0 - eccentricity)).powf(-1.5);
        let mut half_the_excitation_frequency = abs!(particle_norm_spin_vector - mean_motion);
        // If the dynamical tide is excited, compute the planet-dependent stellar dissipation
        // If the dynamical tide is not excited (i.e., equilibrium tide), do nothing since the default value corresponds to the planet scaled dissipation factor
        if half_the_excitation_frequency < particle_norm_spin_vector {
            // The dynamical tide is excited
            if half_the_excitation_frequency < SMOOTHING_FACTOR_DYN_TIDE_COROTATION {
                half_the_excitation_frequency = SMOOTHING_FACTOR_DYN_TIDE_COROTATION;
            }
            let inverse_of_half_the_excitation_frequency = 1. / half_the_excitation_frequency;

            let dynamical = dynamical_scaled_dissipation_factor(
                &params,
                particle.radius,
                particle.tides.parameters.internal.lag_angle,
                inverse_of_half_the_excitation_frequency,
            );
            let host_dependent_scaled_dissipation_factor = match params.tide_composition {
                TideComposition::Dynamical => dynamical,
                // Historical behavior: dissipation from the equilibrium tide
                // is added to the dynamical tide one
                _ => dynamical + equilibrium_scaled_dissipation_factor(&params),
            };

            set_pair_dependent_scaled_dissipation_factor(
                pair_dependent_scaled_dissipation_factor,
                particle.id,
                tidal_host_particle.id,
                host_dependent_scaled_dissipation_factor,
            );
        } else {
            // Equilibrium regime (not dynamical tide).
            // Note: the key must match the one used by the set above —
            // (particle, host) — so that leaving the dynamical regime clears
            // this particle's stale override (and not the host's entry).
            remove_pair_dependent_scaled_dissipation_factor(
                pair_dependent_scaled_dissipation_factor,
                particle.id,
                tidal_host_particle.id,
            );
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
    tide_composition: TideComposition,
    equilibrium_scaled_dissipation_factor: f64,
) -> f64 {
    let key = id * MAX_PARTICLES + depends_on_id;
    match tide_composition {
        TideComposition::Equilibrium => equilibrium_scaled_dissipation_factor,
        TideComposition::Dynamical => {
            if evolution_provides_inverse_tidal_q(evolution) {
                // Outside the excitation regime the dynamical tide does not
                // dissipate at all (no equilibrium fallback in this mode).
                pair_dependent_scaled_dissipation_factor
                    .get(&key)
                    .copied()
                    .unwrap_or(0.)
            } else {
                // The chosen evolution model provides no 1/Q: the dynamical
                // tide cannot be computed and this body does not dissipate.
                0.
            }
        }
        TideComposition::Both => {
            if evolution_provides_inverse_tidal_q(evolution) {
                match pair_dependent_scaled_dissipation_factor.get(&key) {
                    Some(&value) => value,
                    _ => equilibrium_scaled_dissipation_factor, // This happens if it goes out of the dynamical regime into the equilibrium regime
                }
            } else {
                equilibrium_scaled_dissipation_factor
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
/// TIDES
pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    let orthogonal_component_of_the_tidal_force;
    let reference_rscalspin;
    let reference_spin;

    if central_body {
        reference_spin = tidal_host_particle.spin;
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
        reference_spin = particle.spin;
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
    calculate_orthogonal_component_of_the_stellar_tide(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
    calculate_orthogonal_component_of_the_planetary_tide(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
}

/// Orthogonal component of the tide raised IN THE STAR by each companion.
///
/// Gated on the STAR being a CTL central body: the companion only contributes
/// its mass and coordinates, so its own tidal model (or lack of one) is
/// irrelevant. The result is stored per companion.
fn calculate_orthogonal_component_of_the_stellar_tide(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    let host_params = match &tidal_host_particle.tides.effect {
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(params)) => *params,
        _ => return,
    };
    let host_equilibrium_scaled_dissipation_factor =
        equilibrium_scaled_dissipation_factor(&host_params);
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // Defensive: skip degenerate coordinates (no meaningful orbit)
        if particle.tides.parameters.internal.distance == 0. {
            continue;
        }
        // (distance to star)^7
        let distance_7 = particle.tides.parameters.internal.distance.powi(7);

        // - Third line of Equation 5 from Bolmont et al. 2015
        //   This expression has R**10 (instead of R**5 in Eq. 5)
        //   because it uses sigma (i.e., scaled_dissipation_factor)
        //   and not k2$\Delta$t (between k2$\Delta$t and sigma
        //   there is a R**5 factor as shown in Equation 28)
        //   - k2 is love number
        let host_scaled_dissipation_factor = get_pair_dependent_scaled_dissipation_factor_or_else(
            pair_dependent_scaled_dissipation_factor,
            tidal_host_particle.id,
            particle.id,
            tidal_host_particle.evolution,
            host_params.tide_composition,
            host_equilibrium_scaled_dissipation_factor,
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
    }
}

/// Orthogonal component of the tide raised IN EACH PLANET by the star.
///
/// Gated on the PLANET being a CTL orbiting body.
fn calculate_orthogonal_component_of_the_planetary_tide(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        if let TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(params)) =
            &particle.tides.effect
        {
            let params = *params;
            // (distance to star)^7
            let distance_7 = particle.tides.parameters.internal.distance.powi(7);

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
                    params.tide_composition,
                    equilibrium_scaled_dissipation_factor(&params),
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
        }
    }
}

pub fn calculate_radial_component_of_the_tidal_force(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    calculate_radial_component_of_the_stellar_tide(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
    calculate_radial_component_of_the_planetary_tide(
        tidal_host_particle,
        particles,
        more_particles,
        pair_dependent_scaled_dissipation_factor,
    );
}

/// Radial component (conservative + dissipative) of the tide raised IN THE
/// STAR by each companion — the terms of the first line of Equation 5 from
/// Bolmont et al. 2015 that carry the star's radius, love number and sigma.
///
/// Gated on the STAR being a CTL central body; stored per companion.
fn calculate_radial_component_of_the_stellar_tide(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
    pair_dependent_scaled_dissipation_factor: &mut HashMap<usize, f64>,
) {
    let host_params = match &tidal_host_particle.tides.effect {
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(params)) => *params,
        _ => return,
    };
    let host_equilibrium_scaled_dissipation_factor =
        equilibrium_scaled_dissipation_factor(&host_params);
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // Defensive: skip degenerate coordinates (no meaningful orbit)
        if particle.tides.parameters.internal.distance == 0. {
            continue;
        }
        let particle_mass_2 = particle.mass * particle.mass;
        // Conservative part
        let conservative_part = -3.0 * K2 / particle.tides.parameters.internal.distance.powi(7)
            * particle_mass_2
            * tidal_host_particle.radius.powi(5)
            * host_params.love_number;
        // Dissipative part
        let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity
            / particle.tides.parameters.internal.distance.powi(8);
        let host_scaled_dissipation_factor = get_pair_dependent_scaled_dissipation_factor_or_else(
            pair_dependent_scaled_dissipation_factor,
            tidal_host_particle.id,
            particle.id,
            tidal_host_particle.evolution,
            host_params.tide_composition,
            host_equilibrium_scaled_dissipation_factor,
        );
        let term1 =
            particle_mass_2 * tidal_host_particle.radius.powi(10) * host_scaled_dissipation_factor;
        particle
            .tides
            .parameters
            .internal
            .radial_component_of_the_tidal_force_due_to_stellar_tide =
            conservative_part + factor1 * term1;
    }
}

/// Radial component (conservative + dissipative) of the tide raised IN EACH
/// PLANET by the star — the terms of the first line of Equation 5 from
/// Bolmont et al. 2015 that carry the planet's radius, love number and sigma.
///
/// Gated on the PLANET being a CTL orbiting body.
fn calculate_radial_component_of_the_planetary_tide(
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
            let params = *params;
            // Conservative part
            let conservative_part = -3.0 * K2 / particle.tides.parameters.internal.distance.powi(7)
                * host_mass_2
                * particle.radius.powi(5)
                * params.love_number;

            // Dissipative part
            let factor1 = -13.5 * particle.tides.parameters.internal.radial_velocity
                / particle.tides.parameters.internal.distance.powi(8);
            let particle_scaled_dissipation_factor =
                get_pair_dependent_scaled_dissipation_factor_or_else(
                    pair_dependent_scaled_dissipation_factor,
                    particle.id,
                    tidal_host_particle.id,
                    particle.evolution,
                    params.tide_composition,
                    equilibrium_scaled_dissipation_factor(&params),
                );
            let term2 = host_mass_2 * particle.radius.powi(10) * particle_scaled_dissipation_factor;
            // If we consider the star as a point mass (used for denergy_dt calculation):
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass =
                factor1 * term2;
            particle
                .tides
                .parameters
                .internal
                .radial_component_of_the_tidal_force_due_to_planetary_tide = conservative_part
                + particle
                    .tides
                    .parameters
                    .internal
                    .radial_component_of_the_tidal_force_dissipative_part_when_star_as_point_mass;
        }
    }
}

// - Equation 6 from Bolmont et al. 2015, split per tide.
//
// The total pair force is linear in the per-tide radial and orthogonal
// components, so it separates exactly into a stellar-tide part (terms carrying
// the star's radius/love number/sigma and the star's spin) and a planetary-
// tide part (the planet's). Both functions return the force ACTING ON THE
// PLANET; the force on the star is the Newton's-third-law reaction.
fn calculate_tidal_force_for(
    radial_component: f64,
    orthogonal_component: f64,
    spin: &Axes,
    particle: &Particle,
) -> Axes {
    let distance = particle.tides.parameters.internal.distance;
    let position = &particle.tides.coordinates.position;
    let velocity = &particle.tides.coordinates.velocity;
    let factor3 = radial_component
        + orthogonal_component * particle.tides.parameters.internal.radial_velocity / distance;
    let total_tidal_force_x = factor3 * position.x / distance
        + orthogonal_component / distance
            * (spin.y * position.z - spin.z * position.y - velocity.x);
    let total_tidal_force_y = factor3 * position.y / distance
        + orthogonal_component / distance
            * (spin.z * position.x - spin.x * position.z - velocity.y);
    let total_tidal_force_z = factor3 * position.z / distance
        + orthogonal_component / distance
            * (spin.x * position.y - spin.y * position.x - velocity.z);
    Axes::from(
        total_tidal_force_x,
        total_tidal_force_y,
        total_tidal_force_z,
    )
}

/// Force ON THE PLANET due to the tide raised IN THE STAR by this planet
/// (couples to the star's spin).
pub fn calculate_stellar_tidal_force(tidal_host_particle: &Particle, particle: &Particle) -> Axes {
    calculate_tidal_force_for(
        particle
            .tides
            .parameters
            .internal
            .radial_component_of_the_tidal_force_due_to_stellar_tide,
        particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide,
        &tidal_host_particle.spin,
        particle,
    )
}

/// Force ON THE PLANET due to the tide raised IN THE PLANET by the star
/// (couples to the planet's spin).
pub fn calculate_planetary_tidal_force(particle: &Particle) -> Axes {
    calculate_tidal_force_for(
        particle
            .tides
            .parameters
            .internal
            .radial_component_of_the_tidal_force_due_to_planetary_tide,
        particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide,
        &particle.spin,
        particle,
    )
}
