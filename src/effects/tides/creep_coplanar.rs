// Implemented by Gabriel Oliveira Gomes
// Described in Gomes et al (2021): https://ui.adsabs.harvard.edu/abs/2021A%26A...651A..23G/abstract
use super::super::super::Axes;
use super::super::super::Particle;
use super::super::super::constants::{K2, PI, TWO_PI};
use super::super::super::tools;
use super::common::TidalModel;
use super::common::TidesEffect;
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct CreepCoplanarParameters {
    pub uniform_viscosity_coefficient: f64,
}

pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    // Torque expression for studying creep tide tidal despinning (use torque_z = 0 instead to study stat. rotation, makes code run faster)

    let factor = 3.0 / 5.0
        * K2
        * particle.mass
        * tidal_host_particle.mass
        * tidal_host_particle.radius.powi(2);

    let torque_due_to_tides_z = if central_body {
        // Host shape is planet dependent, it needs to be computed for each.
        // For (mainly) co-orbital and circumbinary systems, planet shape cannot be linearly added.
        // The resulting shape is a combination of the shapes caused by each planet, but the
        // way of implementation is not a direct sum of the force components, but actually a
        // composition of the figures of equilibrium caused by each planet.
        // See more discussions in https://arxiv.org/abs/2105.02336
        let consider_tides = true;
        let consider_rotational_flattening = false;
        let tidal_host_shape = calculate_creep_coplanar_shape(
            tidal_host_particle,
            particle,
            consider_tides,
            consider_rotational_flattening,
            central_body,
        );
        factor * tidal_host_shape.y / particle.tides.parameters.internal.distance.powi(3)
    } else {
        // Planet shape is host dependent and it was already computed
        factor * particle.tides.parameters.internal.shape.y
            / particle.tides.parameters.internal.distance.powi(3)
    };

    Axes::from(0.0, 0.0, torque_due_to_tides_z)
}

// Expression for ONLY tides
fn tidal_force_only_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    position: f64,
    sign: f64,
) -> f64 {
    let distance_4 = particle.tides.parameters.internal.distance.powi(4);
    let factor1 = K2 * particle.mass * tidal_host_particle.mass * particle.radius.powi(2);

    particle.tides.parameters.internal.distance
        * (-0.9 * factor1 / distance_4 * particle.tides.parameters.internal.shape.x
            - (3.0 * factor1 / (5.0 * distance_4) * particle.tides.parameters.internal.shape.z))
        + sign
            * (3.0 * factor1 * position
                / (5.0 * distance_4 * particle.tides.parameters.internal.distance)
                * particle.tides.parameters.internal.shape.y)
}

pub fn calculate_tidal_force(tidal_host_particle: &Particle, particle: &Particle) -> Axes {
    // Expression for ONLY tides
    let sign = -1.0;
    let total_tidal_force_x = particle.tides.coordinates.position.x
        / tidal_force_only_tides(
            tidal_host_particle,
            particle,
            particle.tides.coordinates.position.y,
            sign,
        );

    // Expression for ONLY tides
    let sign = 1.0;
    let total_tidal_force_y = particle.tides.coordinates.position.y
        / tidal_force_only_tides(
            tidal_host_particle,
            particle,
            particle.tides.coordinates.position.x,
            sign,
        );

    Axes::from(total_tidal_force_x, total_tidal_force_y, 0.)
}

////////////////////////////////////////////////////////////////////////////////
// Creep coplanar tools
////////////////////////////////////////////////////////////////////////////////
pub fn calculate_creep_coplanar_shapes(
    tidal_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    // It does not compute the host shape since it is planet dependent
    // and it will be computed later when it is needed
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        // particle_uniform_viscosity_coefficient
        let (TidesEffect::CentralBody(TidalModel::CreepCoplanar(params))
        | TidesEffect::OrbitingBody(TidalModel::CreepCoplanar(params))) = &particle.tides.effect
        else {
            continue;
        };
        if params.uniform_viscosity_coefficient == 0.0 {
            continue;
        }

        let consider_tides = true;
        let consider_rotational_flattening = false;
        let central_body = false;
        let shape = calculate_creep_coplanar_shape(
            tidal_host_particle,
            particle,
            consider_tides,
            consider_rotational_flattening,
            central_body,
        );
        particle.tides.parameters.internal.shape = shape;
    }
}

pub fn calculate_creep_coplanar_shape(
    tidal_host_particle: &Particle,
    particle: &Particle,
    consider_tides: bool,
    consider_rotational_flattening: bool,
    central_body: bool,
) -> Axes {
    let gm = tidal_host_particle.mass_g + particle.mass_g;
    let (
        semimajor_axis,
        _perihelion_distance,
        eccentricity,
        _inclination,
        _perihelion_longitude,
        _longitude_of_ascending_node,
        mean_anomaly,
        orbital_period,
    ) = tools::calculate_keplerian_orbital_elements(
        gm,
        particle.tides.coordinates.position,
        particle.tides.coordinates.velocity,
    );

    if central_body {
        let tidal_host_particle_uniform_viscosity_coefficient =
            match &tidal_host_particle.tides.effect {
                TidesEffect::CentralBody(TidalModel::CreepCoplanar(params)) => {
                    params.uniform_viscosity_coefficient
                }
                _ => 0.,
            };
        calculate_particle_shape(
            semimajor_axis,
            eccentricity,
            mean_anomaly,
            orbital_period,
            tidal_host_particle.mass,
            particle.mass,
            tidal_host_particle.radius,
            tidal_host_particle_uniform_viscosity_coefficient,
            tidal_host_particle.spin.z,
            consider_tides,
            consider_rotational_flattening,
        )
    } else {
        let particle_uniform_viscosity_coefficient = match &particle.tides.effect {
            TidesEffect::OrbitingBody(TidalModel::CreepCoplanar(params)) => {
                params.uniform_viscosity_coefficient
            }
            _ => 0.,
        };
        calculate_particle_shape(
            semimajor_axis,
            eccentricity,
            mean_anomaly,
            orbital_period,
            particle.mass,
            tidal_host_particle.mass,
            particle.radius,
            particle_uniform_viscosity_coefficient,
            particle.spin.z,
            consider_tides,
            consider_rotational_flattening,
        )
    }
}

pub fn calculate_particle_shape(
    semimajor_axis: f64,
    eccentricity: f64,
    mean_anomaly: f64,
    orbital_period: f64,
    primary_mass: f64,
    companion_mass: f64,
    radius: f64,
    uniform_viscosity_coefficient: f64,
    primary_spin: f64,
    consider_tides: bool,
    consider_rotational_flattening: bool,
) -> Axes {
    let x_component_of_the_equatorial_prolateness;
    let y_component_of_the_equatorial_prolateness;
    let epsilon_z;

    let relaxation_factor = 3.0 * K2 * primary_mass * primary_mass
        / (8.0 * PI * radius.powi(4) * uniform_viscosity_coefficient);
    let mean_motion = TWO_PI / orbital_period;
    let epsilon_rho_bar = if consider_tides {
        15.0 * companion_mass * radius.powi(3) / (4.0 * primary_mass * semimajor_axis.powi(3))
    } else {
        0.
    };

    let epsilon_z_bar = if consider_rotational_flattening {
        5.0 * primary_spin * primary_spin * radius.powi(3) / (4.0 * K2 * primary_mass)
    } else {
        0.
    };
    // This part is important, adjust properly!

    let pseudo_synchronization_frequency = if relaxation_factor < mean_motion {
        0.9 * mean_motion // 1.0 * mean_motion + 6.0 * eccentricity * eccentricity * mean_motion;
    } else {
        1.0 * mean_motion + 7.0 * eccentricity * eccentricity * mean_motion
        // 1.0 * mean_motion + 8.0 * eccentricity * eccentricity * mean_motion;
    };

    if primary_spin < pseudo_synchronization_frequency {
        // be careful here
        // The Epsilons are taken from Equation 49 of Folonier et al. (2018)
        let companion_mass_to_sum_of_masses_ratio =
            companion_mass / (primary_mass + companion_mass);

        let alpha = 1.0 - (3.0 * companion_mass_to_sum_of_masses_ratio * epsilon_rho_bar);
        let mean_motion_to_relaxation_factor_ratio = mean_motion / relaxation_factor;
        let mean_motion_to_relaxation_factor_ratio_2 =
            mean_motion_to_relaxation_factor_ratio.powi(2);
        // eccentricity_2 = eccentricity.powi(2);
        // alpha_2 = alpha.powi(2);

        let equatorial_prolateness_constant_coefficient = epsilon_rho_bar
            * (1.0 + 1.5 * eccentricity.powi(2)
                - (4.0 * mean_motion_to_relaxation_factor_ratio_2 * eccentricity.powi(2)
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2)));
        let equatorial_prolateness_cosine_coefficient =
            3.0 * epsilon_rho_bar * eccentricity / (1.0 + mean_motion_to_relaxation_factor_ratio_2);
        let equatorial_prolateness_sine_coefficient =
            3.0 * epsilon_rho_bar * eccentricity * mean_motion_to_relaxation_factor_ratio
                / (1.0 + mean_motion_to_relaxation_factor_ratio_2);

        let polar_oblateness_constant_coefficient = epsilon_rho_bar
            * (0.5 + 0.75 * eccentricity.powi(2))
            + epsilon_z_bar
                * (1.0
                    + 12.0
                        * eccentricity.powi(2)
                        * (1.0 + alpha * mean_motion_to_relaxation_factor_ratio_2)
                        / ((1.0 + mean_motion_to_relaxation_factor_ratio_2)
                            * (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2))
                    + 2.0
                        * (1.0 - alpha)
                        * (1.0 - alpha)
                        * mean_motion_to_relaxation_factor_ratio_2
                        * eccentricity.powi(2)
                        / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2));
        let polar_oblateness_cosine_coefficient = (1.5 * epsilon_rho_bar * eccentricity
            / (1.0 + mean_motion_to_relaxation_factor_ratio_2))
            * (1.0
                - (16.0
                    * companion_mass_to_sum_of_masses_ratio
                    * mean_motion_to_relaxation_factor_ratio_2
                    * epsilon_z_bar
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2)));
        let polar_oblateness_sine_coefficient = (1.5 * epsilon_rho_bar * eccentricity
            / (1.0 + mean_motion_to_relaxation_factor_ratio_2))
            * (mean_motion_to_relaxation_factor_ratio
                + (8.0
                    * companion_mass_to_sum_of_masses_ratio
                    * mean_motion_to_relaxation_factor_ratio
                    * (1.0 - alpha * mean_motion_to_relaxation_factor_ratio_2)
                    * epsilon_z_bar
                    / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2)));

        let delta_constant_coefficient =
            (3.0 * mean_motion_to_relaxation_factor_ratio * eccentricity.powi(2)
                / (1.0 + mean_motion_to_relaxation_factor_ratio_2))
                * (2.0 + (1.0 + alpha) * mean_motion_to_relaxation_factor_ratio_2)
                / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2);
        let delta_cosine_coefficient = -2.0 * mean_motion_to_relaxation_factor_ratio * eccentricity
            / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2);
        let delta_sine_coefficient =
            -2.0 * eccentricity * mean_motion_to_relaxation_factor_ratio_2 * alpha
                / (1.0 + alpha * alpha * mean_motion_to_relaxation_factor_ratio_2);

        let ep_rho = equatorial_prolateness_constant_coefficient
            + equatorial_prolateness_cosine_coefficient * cos!(mean_anomaly)
            + equatorial_prolateness_sine_coefficient * sin!(mean_anomaly);

        epsilon_z = polar_oblateness_constant_coefficient
            + polar_oblateness_cosine_coefficient * cos!(mean_anomaly)
            + polar_oblateness_sine_coefficient * sin!(mean_anomaly);

        let delta = delta_constant_coefficient
            + delta_cosine_coefficient * cos!(mean_anomaly)
            + delta_sine_coefficient * sin!(mean_anomaly);

        let delta2 = 2.0 * delta;

        x_component_of_the_equatorial_prolateness = ep_rho * cos!(delta2);
        y_component_of_the_equatorial_prolateness = ep_rho * sin!(delta2);
    } else {
        let cayley = calculate_cayley_coefficients(eccentricity);
        let mut true_anomaly = calculate_true_anomaly(eccentricity, mean_anomaly);
        if true_anomaly < 0.0 {
            true_anomaly += 2.0 * PI;
        }

        let mut sum_x_component_of_the_equatorial_prolateness;
        let mut sum_y_component_of_the_equatorial_prolateness;
        let mut sum_epsilon_z;
        let mut equatorial_angles;
        let mut polar_angles;

        sum_x_component_of_the_equatorial_prolateness = 0.0;
        sum_y_component_of_the_equatorial_prolateness = 0.0;
        sum_epsilon_z = 0.0;

        // for (k,factor) in cayley[1].iter().enumerate()

        for j in 0..15 {
            let k = j as f64;
            equatorial_angles = (9.0 - k) * mean_anomaly - 2.0 * true_anomaly;
            polar_angles = (-7.0 + k) * mean_anomaly;

            // common block (relaxation_factor * cayley[1][j] / (relaxation_factor.powi(2) + (2.0*primary_spin-(9.0-k)*mean_motion).powi(2)))

            sum_x_component_of_the_equatorial_prolateness += relaxation_factor * cayley[1][j]
                / (relaxation_factor.powi(2)
                    + (2.0 * primary_spin - (9.0 - k) * mean_motion).powi(2))
                * (relaxation_factor * cos!(equatorial_angles)
                    - (2.0 * primary_spin - (9.0 - k) * mean_motion) * sin!(equatorial_angles));
            sum_y_component_of_the_equatorial_prolateness += relaxation_factor * cayley[1][j]
                / (relaxation_factor.powi(2)
                    + (2.0 * primary_spin - (9.0 - k) * mean_motion).powi(2))
                * (relaxation_factor * sin!(equatorial_angles)
                    + (2.0 * primary_spin - (9.0 - k) * mean_motion) * cos!(equatorial_angles));
            sum_epsilon_z += cayley[0][j]
                / ((k - 7.0).powi(2) * mean_motion.powi(2) + relaxation_factor.powi(2))
                * ((k - 7.0).powi(2) * mean_motion.powi(2) * sin!(polar_angles)
                    + relaxation_factor.powi(2) * cos!(polar_angles));
            //sum_x_component_of_the_equatorial_prolateness += 0.0;
            //sum_y_component_of_the_equatorial_prolateness += 0.0;
            //sum_epsilon_z += 0.0;
        }

        // println!("{}", 0.5 * sum_y_component_of_the_equatorial_prolateness / sum_x_component_of_the_equatorial_prolateness * 180.0 / PI);

        sum_epsilon_z = (sum_epsilon_z * epsilon_rho_bar / 2.0) + epsilon_z_bar; // In case tides + rotation are taken into account

        //sum_epsilon_z = epsilon_z_bar; // In case only rotational flattenings are taken into account

        x_component_of_the_equatorial_prolateness =
            sum_x_component_of_the_equatorial_prolateness * epsilon_rho_bar;
        y_component_of_the_equatorial_prolateness =
            sum_y_component_of_the_equatorial_prolateness * epsilon_rho_bar;
        epsilon_z = sum_epsilon_z;
    }

    Axes::from(
        x_component_of_the_equatorial_prolateness,
        y_component_of_the_equatorial_prolateness,
        epsilon_z,
    )
}

fn calculate_cayley_coefficients(e: f64) -> [[f64; 15]; 2] {
    let mut cayley = [[0.; 15]; 2];

    let e2 = e.powi(2);
    let e3 = e * e2;
    let e4 = e2 * e2;
    let e5 = e4 * e;
    let e6 = e5 * e;
    let e7 = e6 * e;

    // Coefficients E_(0,k) of FM2015 arvix

    cayley[0][0] = 432_091.0 / 30_720.0 * e7;
    cayley[0][1] = 3167.0 / 320.0 * e6;
    cayley[0][2] = 1773.0 / 256.0 * e5 - 4987.0 / 6144.0 * e7;
    cayley[0][3] = 77.0 / 16.0 * e4 + 129.0 / 160.0 * e6;
    cayley[0][4] = 53.0 / 16.0 * e3 + 393.0 / 256.0 * e5 + 24_753.0 / 10_240.0 * e7;
    cayley[0][5] = 9.0 / 4.0 * e2 + 7.0 / 4.0 * e4 + 141.0 / 64.0 * e6;
    cayley[0][6] = 3.0 / 2.0 * e + 27.0 / 16.0 * e3 - 261.0 / 128.0 * e5 + 14_309.0 / 6144.0 * e7;
    cayley[0][7] = 1.0 + 3.0 / 2.0 * e2 + 15.0 / 8.0 * e4 - 35.0 / 16.0 * e6;
    cayley[0][8] = cayley[0][6];
    cayley[0][9] = cayley[0][5];
    cayley[0][10] = cayley[0][4];
    cayley[0][11] = cayley[0][3];
    cayley[0][12] = cayley[0][2];
    cayley[0][13] = cayley[0][1];
    cayley[0][14] = cayley[0][0];

    // Coefficients E_(2,k) of FM2015 arvix

    cayley[1][0] = 12_144_273.0 / 71_680.0 * e7;
    cayley[1][1] = 73_369.0 / 720.0 * e6;
    cayley[1][2] = 228_347.0 / 3840.0 * e5 - 3_071_075.0 / 18_432.0 * e7;
    cayley[1][3] = 533.0 / 16.0 * e4 - 13_827.0 / 160.0 * e6;
    cayley[1][4] = 845.0 / 48.0 * e3 - 32_525.0 / 768.0 * e5 + 208_225.0 / 6144.0 * e7;
    cayley[1][5] = 17.0 / 2.0 * e2 - 115.0 / 6.0 * e4 + 601.0 / 48.0 * e6;
    cayley[1][6] = 7.0 / 2.0 * e - 123.0 / 16.0 * e3 + 489.0 / 128.0 * e5 - 1763.0 / 2048.0 * e7;
    cayley[1][7] = 1.0 - 5.0 / 2.0 * e2 + 13.0 / 16.0 * e4 - 35.0 / 288.0 * e6;
    cayley[1][8] = -1.0 / 2.0 * e + 1.0 / 16.0 * e3 - 5.0 / 384.0 * e5 - 143.0 / 18_432.0 * e7;
    cayley[1][9] = 0.0;
    cayley[1][10] = 1.0 / 48.0 * e3 + 11.0 / 768.0 * e5 + 313.0 / 30_720.0 * e7;
    cayley[1][11] = 1.0 / 24.0 * e4 + 7.0 / 240.0 * e6;
    cayley[1][12] = 81.0 / 1280.0 * e5 + 81.0 / 2048.0 * e7;
    cayley[1][13] = 4.0 / 45.0 * e6;
    cayley[1][14] = 15_625.0 / 129_024.0 * e7;

    cayley
}

fn calculate_true_anomaly(eccentricity: f64, mean_anomaly: f64) -> f64 {
    let mut difference = 10.0;
    let precision = 1.0e-13;
    let mut eccentric_anomaly_0 = mean_anomaly;

    while difference > precision {
        let eccentric_anomaly_1 = mean_anomaly + eccentricity * sin!(eccentric_anomaly_0);
        difference = abs!(eccentric_anomaly_1 - eccentric_anomaly_0);
        eccentric_anomaly_0 = eccentric_anomaly_1;
    }

    2.0 * atan!(
        ((1.0 + eccentricity) / (1.0 - eccentricity)).powf(1.0 / 2.0)
            * tan!(eccentric_anomaly_0 / 2.0)
    )
}
