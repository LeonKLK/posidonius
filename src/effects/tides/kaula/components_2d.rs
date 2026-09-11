use super::{
    KaulaParameters, Orbit, alpha_pqkj, calculate_2d_constant, calculate_base_constant,
    select_eccentricty_order_q,
};

use crate::Particle;
use crate::tools::KeplerianElements;


/// 2D (coplanar) kaula force components on `tidal_deformed_body` due to `tidal_perturber`.
/// Roles: see `kaula.rs`. The orbit is always the planet's.
pub fn calculate_2d_tidal_force_components(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    orbit: &Orbit,
    central_body: bool,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let kaula = tidal_deformed_body.tides.get_kaula();

    let orthogonal_constant = orthogonal_constant(
        tidal_deformed_body,
        tidal_perturber,
        orbit,
        keplerian_elements.semi_major_axis,
        central_body,
    );
    let radial_constant = radial_constant(
        tidal_deformed_body,
        tidal_perturber,
        orbit,
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
// Bug exist for now: imk2_2200 should be rek2_2200
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

/// Heliocentric distance entering the 2D prefactors. For the stellar tide the sign is flipped
/// so that the force projected with the planet's angles points the right way on the star.
fn signed_heliocentric_distance(orbit: &Orbit, central_body: bool) -> f64 {
    if central_body {
        -orbit.heliocentric_distance
    } else {
        orbit.heliocentric_distance
    }
}

fn orthogonal_constant(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    orbit: &Orbit,
    semi_major_axis: f64,
    central_body: bool,
) -> f64 {
    calculate_2d_constant(
        tidal_deformed_body,
        tidal_perturber,
        orbit,
        signed_heliocentric_distance(orbit, central_body),
        semi_major_axis,
    )
}

fn radial_constant(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    orbit: &Orbit,
    semi_major_axis: f64,
    central_body: bool,
) -> f64 {
    -calculate_base_constant(
        tidal_deformed_body,
        tidal_perturber,
        signed_heliocentric_distance(orbit, central_body),
        semi_major_axis,
    )
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
    let n_q = q_max - q_min;

    // With p = k = 0 the phase alpha_pqkj is (j - q) * mean_anomaly (the argument of perihelion
    // term vanishes; kept in the call so the angle is formed as in `alpha_pqkj`). The double sum
    // over (q, j) then separates with the angle-addition formulas:
    //   sum_j g_j cos((j - q) M) = cos(q M) sum_j g_j cos(j M) + sin(q M) sum_j g_j sin(j M)
    //   sum_j g_j sin((j - q) M) = cos(q M) sum_j g_j sin(j M) - sin(q M) sum_j g_j cos(j M)
    // so the four q-independent sums are formed once and each q costs O(1) instead of O(n_q).
    let g_20 = &kaula.polynomials.eccentricity_function_g_2pq[0][q_min..q_max];
    let g_21 = &kaula.polynomials.eccentricity_function_g_2pq[1][q_min..q_max];
    let mut cos_t = [0.0_f64; 15];
    let mut sin_t = [0.0_f64; 15];
    let (mut a_20, mut b_20, mut a_21, mut b_21) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for t in 0..n_q {
        let alpha = alpha_pqkj(0, 0, 0, t, mean_anomaly, argument_perihelion); // t * M
        cos_t[t] = cos!(alpha);
        sin_t[t] = sin!(alpha);
        a_20 += g_20[t] * sin_t[t];
        b_20 += g_20[t] * cos_t[t];
        a_21 += g_21[t] * sin_t[t];
        b_21 += g_21[t] * cos_t[t];
    }

    // For q in the range defined by the eccentricity
    for q in 0..n_q {
        let (g_20q, g_21q) = (g_20[q], g_21[q]);
        let tmp_q = q_min + q;
        let rek2_201q = kaula.love_numbers.real(0, 1, tmp_q);
        let imk2_201q = kaula.love_numbers.imaginary(0, 1, tmp_q);
        let rek2_220q = kaula.love_numbers.real(2, 0, tmp_q);
        let imk2_220q = kaula.love_numbers.imaginary(2, 0, tmp_q);

        let (cos_q, sin_q) = (cos_t[q], sin_t[q]);
        // sum_j g_2pj cos((j - q) M) and sum_j g_2pj sin((j - q) M) for p = 0 and p = 1
        let c_20 = cos_q * b_20 + sin_q * a_20;
        let s_20 = cos_q * a_20 - sin_q * b_20;
        let c_21 = cos_q * b_21 + sin_q * a_21;
        let s_21 = cos_q * a_21 - sin_q * b_21;

        let sum_over_j_1 = c_21 * rek2_201q - s_21 * imk2_201q;
        let sum_over_j_2 = s_20 * rek2_220q + c_20 * imk2_220q;
        let sum_over_j_3 = c_20 * rek2_220q - s_20 * imk2_220q;
        // secular part: only the j = q terms (alpha = 0)
        let sum_over_j_1_secular = g_21q * rek2_201q;
        let sum_over_j_2_secular = g_20q * imk2_220q;
        let sum_over_j_3_secular = g_20q * rek2_220q;

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
