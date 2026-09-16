use super::{
    KaulaParameters, alpha_pqkj, calculate_2d_constant, calculate_base_constant,
    select_eccentricty_order_q,
};

use crate::Particle;
use crate::tools::KeplerianElements;

use itertools::izip;

pub fn calculate_2d_tidal_force_components(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let kaula = particle.tides.get_kaula();

    let orthogonal_constant = orthogonal_constant(
        tidal_host_particle,
        particle,
        keplerian_elements.semi_major_axis,
        central_body,
    );
    let radial_constant = radial_constant(
        tidal_host_particle,
        particle,
        keplerian_elements.semi_major_axis,
        central_body,
    );

    let (normal_component, normal_component_secular) = (0., 0.);

    let (
        (orthogonal_component, orthogonal_component_secular),
        (radial_component, radial_component_secular),
    ) = if keplerian_elements.eccentricity == 0.0 {
        // if the inclination between the equatorial plane and the orbital plane is negligible
        zero_eccentricity_components(kaula)
    } else {
        calculate_2d_components(kaula, keplerian_elements)
    };

    (
        (
            normal_component,
            orthogonal_component * orthogonal_constant,
            radial_component * radial_constant,
        ),
        (
            normal_component_secular,
            orthogonal_component_secular * orthogonal_constant,
            radial_component_secular * radial_constant,
        ),
    )
}

// TODO document this function
fn zero_eccentricity_components(kaula: &KaulaParameters) -> ((f64, f64), (f64, f64)) {
    // If circular coplanar orbit
    let rek2_2010 = kaula.love_numbers.real(0, 1, 0);
    let imk2_2200 = kaula.love_numbers.imaginary(2, 0, 0);

    let orthogonal_force = (3. / 2.) * imk2_2200;
    let radial_force = (3. / 4.) * rek2_2010 + (9. / 4.) * imk2_2200;

    // forces and secular forces are the same
    (
        (orthogonal_force, orthogonal_force),
        (radial_force, radial_force),
    )
}

fn orthogonal_constant(
    tidal_host_particle: &Particle,
    particle: &Particle,
    semi_major_axis: f64,
    central_body: bool,
) -> f64 {
    if central_body {
        // Stellar tide (called from tides.rs with the arguments swapped):
        //   `tidal_host_particle` is the PLANET (perturber, orbiting body),
        //   `particle` is the STAR (tidally deformed body).
        // The prefactor needs the star's radius and the planet's mass; the angle sin(theta)
        // and the heliocentric distance are those of the planet on its orbit.
        calculate_2d_constant(
            particle,            // deformed body: star
            tidal_host_particle, // perturber: planet
            tidal_host_particle.tides.coordinates.position,
            -tidal_host_particle.heliocentric_distance,
            semi_major_axis,
        )
    } else {
        // Planetary tide: `particle` is the PLANET (deformed body), `tidal_host_particle` the STAR.
        calculate_2d_constant(
            particle,            // deformed body: planet
            tidal_host_particle, // perturber: star
            particle.tides.coordinates.position,
            particle.heliocentric_distance,
            semi_major_axis,
        )
    }
}

fn radial_constant(
    tidal_host_particle: &Particle,
    particle: &Particle,
    semi_major_axis: f64,
    central_body: bool,
) -> f64 {
    if central_body {
        // Stellar tide: `particle` is the STAR (deformed body), `tidal_host_particle` the PLANET.
        -calculate_base_constant(
            particle,            // deformed body: star
            tidal_host_particle, // perturber: planet
            -tidal_host_particle.heliocentric_distance,
            semi_major_axis,
        )
    } else {
        // Planetary tide: `particle` is the PLANET (deformed body), `tidal_host_particle` the STAR.
        -calculate_base_constant(
            particle,            // deformed body: planet
            tidal_host_particle, // perturber: star
            particle.heliocentric_distance,
            semi_major_axis,
        )
    }
}

fn calculate_2d_components(
    kaula: &KaulaParameters,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64), (f64, f64)) {
    // keplerian elements
    let (
        _semi_major_axis,
        _perihelion_distance,
        eccentricity,
        _inclination,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        _orbital_period,
    ) = keplerian_elements.unpack();

    let argument_perihelion = longitude_perihelion - longitude_of_ascending_node;
    let mut orthogonal_sum_over_q = 0.;
    let mut orthogonal_sum_over_q_secular = 0.;
    let mut radial_sum_over_q = 0.;
    let mut radial_sum_over_q_secular = 0.;
    let (q_min, q_max) = select_eccentricty_order_q(eccentricity);

    // For q in the range defined by the eccentricity
    for (q, (g_20q, g_21q)) in izip!(
        kaula.polynomials.eccentricity_function_g_2pq[0],
        kaula.polynomials.eccentricity_function_g_2pq[1],
    )
    .take(q_max)
    .skip(q_min)
    .enumerate()
    {
        let tmp_q = q_min + q;
        let rek2_201q = kaula.love_numbers.real(0, 1, tmp_q);
        let imk2_201q = kaula.love_numbers.imaginary(0, 1, tmp_q);
        let rek2_220q = kaula.love_numbers.real(2, 0, tmp_q);
        let imk2_220q = kaula.love_numbers.imaginary(2, 0, tmp_q);

        let mut sum_over_j_1 = 0.;
        let mut sum_over_j_2 = 0.;
        let mut sum_over_j_3 = 0.;

        let mut sum_over_j_1_secular = 0.;
        let mut sum_over_j_2_secular = 0.;
        let mut sum_over_j_3_secular = 0.;

        for (j, (g_20j, g_21j)) in izip!(
            kaula.polynomials.eccentricity_function_g_2pq[0],
            kaula.polynomials.eccentricity_function_g_2pq[1]
        )
        .take(q_max)
        .skip(q_min)
        .enumerate()
        {
            let alpha_qj = alpha_pqkj(
                0,
                q_min + q,
                0,
                q_min + j,
                mean_anomaly,
                argument_perihelion,
            );

            sum_over_j_1 += g_21j * (cos!(alpha_qj) * rek2_201q - sin!(alpha_qj) * imk2_201q);
            sum_over_j_2 += g_20j * (sin!(alpha_qj) * rek2_220q + cos!(alpha_qj) * imk2_220q);
            sum_over_j_3 += g_20j * (cos!(alpha_qj) * rek2_220q - sin!(alpha_qj) * imk2_220q);
            if q == j {
                sum_over_j_1_secular += g_21j * rek2_201q;
                sum_over_j_2_secular += g_20j * imk2_220q;
                sum_over_j_3_secular += g_20j * rek2_220q;
            }
        }

        orthogonal_sum_over_q += g_20q * (3. / 2.) * sum_over_j_2;
        orthogonal_sum_over_q_secular += g_20q * (3. / 2.) * sum_over_j_2_secular;

        radial_sum_over_q += (3. / 4.) * g_21q * sum_over_j_1 + (9. / 4.) * g_20q * sum_over_j_3;
        radial_sum_over_q_secular +=
            (3. / 4.) * g_21q * sum_over_j_1_secular + (9. / 4.) * g_20q * sum_over_j_3_secular;
    }

    (
        (orthogonal_sum_over_q, orthogonal_sum_over_q_secular),
        (radial_sum_over_q, radial_sum_over_q_secular),
    )
}
