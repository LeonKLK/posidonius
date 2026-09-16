use crate::Particle;
use crate::tools::KeplerianElements;

use super::{alpha_pqkj, calculate_base_constant, select_eccentricty_order_q};
use itertools::izip;

pub fn calculate_3d_tidal_force_components(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    let (orthogonal_component, orthogonal_component_secular) =
        calculate_orthogonal_component_of_the_tidal_force_3d(
            tidal_host_particle,
            particle,
            keplerian_elements,
            central_body,
        );
    let (radial_component, radial_component_secular) =
        calculate_radial_component_of_the_tidal_force_3d(
            tidal_host_particle,
            particle,
            keplerian_elements,
            central_body,
        );
    let (normal_component, normal_component_secular) =
        calculate_normal_component_of_the_tidal_force_3d(
            tidal_host_particle,
            particle,
            keplerian_elements,
            central_body,
        );

    (
        (normal_component, orthogonal_component, radial_component),
        (
            normal_component_secular,
            orthogonal_component_secular,
            radial_component_secular,
        ),
    )
}

fn compute_phase_alpha_1(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    spin_angle: f64,
    argument_perihelion: f64,
    longitude_ascending_node: f64,
    heliocentric_varphi: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(2 * k - 2 * p + q - j) * mean_anomaly
        + f64!(2 * (k - p)) * argument_perihelion
        + spin_angle
        - longitude_ascending_node
        - heliocentric_varphi
}

fn compute_phase_alpha_2(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    spin_angle: f64,
    argument_perihelion: f64,
    longitude_ascending_node: f64,
    heliocentric_varphi: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(2 * k - 2 * p + q - j) * mean_anomaly + f64!(2 * (k - p)) * argument_perihelion
        - spin_angle
        + longitude_ascending_node
        + heliocentric_varphi
}

// Calculate the phases of the 2-Kaula transformed tidal forces
fn compute_phase_beta(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    spin_angle: f64,
    argument_perihelion: f64,
    longitude_ascending_node: f64,
    heliocentric_varphi: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(4 - 2 * p - 2 * k + q - j) * mean_anomaly
        + f64!(4 - 2 * (k + p)) * argument_perihelion
        + spin_angle
        - longitude_ascending_node
        + heliocentric_varphi
}

fn calculate_normal_component_of_the_tidal_force_3d(
    tidal_host_particle: &Particle,
    particle: &Particle,
    keplerian_elements: &KeplerianElements,
    central_body: bool,
) -> (f64, f64) {
    let kaula = particle.tides.get_kaula();

    // Keplerian elements
    let (
        semi_major_axis,
        _perihelion_distance,
        eccentricity,
        _inclination,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        _orbital_period,
    ) = keplerian_elements.unpack();
    let argument_perihelion = longitude_perihelion - longitude_of_ascending_node;

    // 3D case
    let heliocentric_varphi = 0.;
    let spin_angle = 0.;

    let (q_min, q_max) = select_eccentricty_order_q(eccentricity);

    let mut sum_over_p = 0.;
    let mut sum_over_p_s = 0.;
    for (p, (f_20p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[0],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let mut sum_over_q = 0.;
        let mut sum_over_q_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let tmp_q = q_min + q;
            let rek2_20pq = kaula.love_numbers.real(0, p, tmp_q);
            let imk2_20pq = kaula.love_numbers.imaginary(0, p, tmp_q);

            let mut sum_over_k = 0.;
            let mut sum_over_k_s = 0.;

            for (k, (f_21k, sum_g_2kj)) in izip!(
                kaula.polynomials.inclination_function_f_2mp[1],
                kaula.polynomials.eccentricity_function_g_2pq
            )
            .enumerate()
            {
                let mut sum_over_j = 0.;
                let mut sum_over_j_s = 0.;

                for (j, g_2kj) in sum_g_2kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_beta = compute_phase_beta(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle, // 0.0
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi, // 0.0
                    );
                    let phase_alpha_1 = compute_phase_alpha_1(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle, // 0.0
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi, // 0.0
                    );

                    let cos_alpha_1 = cos!(phase_alpha_1);
                    let sin_alpha_1 = sin!(phase_alpha_1);
                    let cos_beta = cos!(phase_beta);
                    let sin_beta = sin!(phase_beta);

                    sum_over_j += g_2kj
                        * (0.5 * (cos_alpha_1 * rek2_20pq - sin_alpha_1 * imk2_20pq)
                            + 1.5 * (cos_beta * rek2_20pq - sin_beta * imk2_20pq));
                    sum_over_j_s += g_2kj * (0.5 * rek2_20pq + 1.5 * rek2_20pq);
                }
                sum_over_k += f_21k * sum_over_j;
                sum_over_k_s += f_21k * sum_over_j_s;
            }
            sum_over_q += g_2pq * sum_over_k;
            sum_over_q_s += g_2pq * sum_over_k_s;
        }
        sum_over_p += f_20p * sum_over_q;
        sum_over_p_s += f_20p * sum_over_q_s;
    }
    let term_m0 = sum_over_p;
    let term_m0_s = sum_over_p_s;

    sum_over_p = 0.;
    sum_over_p_s = 0.;
    for (p, (f_21p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[1],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let mut sum_over_q = 0.;
        let mut sum_over_q_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let tmp_q = q_min + q;
            let rek2_21pq = kaula.love_numbers.real(1, p, tmp_q);
            let imk2_21pq = kaula.love_numbers.imaginary(1, p, tmp_q);

            let mut sum_over_k_m2 = 0.;
            let mut sum_over_k_m0 = 0.;
            let mut sum_over_k_m2_s = 0.;
            let mut sum_over_k_m0_s = 0.;

            for (k, (f_20k, f_22k, sum_g_2kj)) in izip!(
                kaula.polynomials.inclination_function_f_2mp[0],
                kaula.polynomials.inclination_function_f_2mp[2],
                kaula.polynomials.eccentricity_function_g_2pq
            )
            .enumerate()
            {
                // sum over p for each m
                let mut sum_over_j_p2 = 0.;
                let mut sum_over_j_p0 = 0.;
                let mut sum_over_j_p2_s = 0.;
                let mut sum_over_j_p0_s = 0.;

                for (j, g_2kj) in sum_g_2kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_alpha_1 = compute_phase_alpha_1(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );
                    let phase_alpha_2 = compute_phase_alpha_2(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );

                    let cos_alpha_1 = cos!(phase_alpha_1);
                    let sin_alpha_1 = sin!(phase_alpha_1);
                    let cos_alpha_2 = cos!(phase_alpha_2);
                    let sin_alpha_2 = sin!(phase_alpha_2);

                    sum_over_j_p2 += g_2kj * (cos_alpha_1 * rek2_21pq - sin_alpha_1 * imk2_21pq);
                    sum_over_j_p0 += g_2kj * (cos_alpha_2 * rek2_21pq - sin_alpha_2 * imk2_21pq);
                    sum_over_j_p2_s += g_2kj * rek2_21pq;
                    sum_over_j_p0_s += g_2kj * rek2_21pq;
                }
                sum_over_k_m2 += f_22k * sum_over_j_p2;
                sum_over_k_m0 += f_20k * sum_over_j_p0;
                sum_over_k_m2_s += f_22k * sum_over_j_p2_s;
                sum_over_k_m0_s += f_20k * sum_over_j_p0_s;
            }
            sum_over_q += g_2pq * (0.5 * sum_over_k_m2 - 3. * sum_over_k_m0);
            sum_over_q_s += g_2pq * (0.5 * sum_over_k_m2_s - 3. * sum_over_k_m0_s);
        }
        sum_over_p += f_21p * sum_over_q;
        sum_over_p_s += f_21p * sum_over_q_s;
    }
    let term_m1 = (1. / 3.) * sum_over_p;
    let term_m1_s = (1. / 3.) * sum_over_p_s;

    sum_over_p = 0.;
    sum_over_p_s = 0.;
    for (p, (f_22p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[2],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let mut sum_over_q = 0.;
        let mut sum_over_q_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let tmp_q = q_min + q;
            let rek2_22pq = kaula.love_numbers.real(2, p, tmp_q);
            let imk2_22pq = kaula.love_numbers.imaginary(2, p, tmp_q);

            let mut sum_over_k = 0.;
            let mut sum_over_k_s = 0.;

            for (k, (f_21k, sum_g_2kj)) in izip!(
                kaula.polynomials.inclination_function_f_2mp[1],
                kaula.polynomials.eccentricity_function_g_2pq
            )
            .enumerate()
            {
                let mut sum_over_j = 0.;
                let mut sum_over_j_s = 0.;

                for (j, g_2kj) in sum_g_2kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_alpha_2 = compute_phase_alpha_2(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );

                    let cos_alpha_2 = cos!(phase_alpha_2);
                    let sin_alpha_2 = sin!(phase_alpha_2);

                    sum_over_j += g_2kj * (cos_alpha_2 * rek2_22pq - sin_alpha_2 * imk2_22pq);
                    sum_over_j_s += g_2kj * rek2_22pq;
                }
                sum_over_k += f_21k * sum_over_j;
                sum_over_k_s += f_21k * sum_over_j_s;
            }
            sum_over_q += g_2pq * sum_over_k;
            sum_over_q_s += g_2pq * sum_over_k_s;
        }
        sum_over_p += f_22p * sum_over_q;
        sum_over_p_s += f_22p * sum_over_q_s;
    }
    let term_m2 = -(1. / 6.) * sum_over_p;
    let term_m2_s = -(1. / 6.) * sum_over_p_s;

    let cste = if central_body {
        // Stellar tide: `particle` is the STAR (deformed body, radius^5) and
        // `tidal_host_particle` the PLANET (perturber, mass^2, heliocentric distance).
        calculate_base_constant(
            particle,            // deformed body: star
            tidal_host_particle, // perturber: planet
            tidal_host_particle.heliocentric_distance,
            semi_major_axis,
        )
    } else {
        // Planetary tide: `particle` is the PLANET (deformed body), `tidal_host_particle` the STAR.
        calculate_base_constant(
            particle,
            tidal_host_particle,
            particle.heliocentric_distance,
            semi_major_axis,
        )
    };

    let normal_force = cste * (term_m0 + term_m1 + term_m2);
    let normal_force_secular = cste * (term_m0_s + term_m1_s + term_m2_s);

    (normal_force, normal_force_secular)
}

fn compute_phase_alpha_3(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    spin_angle: f64,
    argument_perihelion: f64,
    longitude_ascending_node: f64,
    heliocentric_varphi: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(2 * k - 2 * p + q - j - 1) * mean_anomaly
        + f64!(2 * (k - p - 1)) * argument_perihelion
        + spin_angle
        - longitude_ascending_node
        - heliocentric_varphi
}

fn compute_phase_alpha_4(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    spin_angle: f64,
    argument_perihelion: f64,
    longitude_ascending_node: f64,
    heliocentric_varphi: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(2 * k - 2 * p + q - j - 1) * mean_anomaly + f64!(2 * (k - p - 1)) * argument_perihelion
        - spin_angle
        + longitude_ascending_node
        - heliocentric_varphi
}

// Ortho-radial (the e_{\varphi}) component of tidal force
fn calculate_orthogonal_component_of_the_tidal_force_3d(
    tidal_host_particle: &Particle,
    particle: &Particle,
    keplerian_elements: &KeplerianElements,
    central_body: bool,
) -> (f64, f64) {
    let kaula = particle.tides.get_kaula();
    // Keplerian elements
    let (
        semi_major_axis,
        _perihelion_distance,
        eccentricity,
        _inclination,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        _orbital_period,
    ) = keplerian_elements.unpack();

    let argument_perihelion = longitude_perihelion - longitude_of_ascending_node;
    // Cartesian angles
    // summation
    let heliocentric_varphi = 0.;
    let spin_angle = 0.;

    let (q_min, q_max) = select_eccentricty_order_q(eccentricity);
    let mut sum_over_p;
    let mut sum_over_p_s;

    sum_over_p = 0.;
    sum_over_p_s = 0.;
    let c1 = 1. / 6.;
    for (p, (f_21p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[1],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let mut sum_over_q = 0.;
        let mut sum_over_q_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let tmp_q = q_min + q;
            let rek2_21pq = kaula.love_numbers.real(1, p, tmp_q);
            let imk2_21pq = kaula.love_numbers.imaginary(1, p, tmp_q);

            let mut sum_over_kterm1 = 0.;
            let mut sum_over_kterm2 = 0.;
            let mut sum_over_kterm1_s = 0.;
            let mut sum_over_kterm2_s = 0.;

            for (k, (f_30k, f_32k, sum_g_3kj)) in izip!(
                kaula.polynomials.inclination_function_f_3mp[0],
                kaula.polynomials.inclination_function_f_3mp[2],
                kaula.polynomials.eccentricity_function_g_3pq
            )
            .enumerate()
            {
                // Sum over p for each m
                let mut sum_over_jterm1 = 0.;
                let mut sum_over_jterm2 = 0.;
                let mut sum_over_jterm1_s = 0.;
                let mut sum_over_jterm2_s = 0.;

                for (j, g_3kj) in sum_g_3kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_alpha_3 = compute_phase_alpha_3(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );
                    let phase_alpha_4 = compute_phase_alpha_4(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );

                    let cos_alpha_1 = cos!(phase_alpha_3);
                    let sin_alpha_1 = sin!(phase_alpha_3);
                    let cos_alpha_2 = cos!(phase_alpha_4);
                    let sin_alpha_2 = sin!(phase_alpha_4);

                    sum_over_jterm1 += g_3kj * (sin_alpha_1 * rek2_21pq - cos_alpha_1 * imk2_21pq);
                    sum_over_jterm2 += g_3kj * (sin_alpha_2 * rek2_21pq - cos_alpha_2 * imk2_21pq);
                    sum_over_jterm1_s -= g_3kj * imk2_21pq;
                    sum_over_jterm2_s -= g_3kj * imk2_21pq;
                }
                sum_over_kterm1 += f_32k * sum_over_jterm1;
                sum_over_kterm2 += f_30k * sum_over_jterm2;
                sum_over_kterm1_s += f_32k * sum_over_jterm1_s;
                sum_over_kterm2_s += f_30k * sum_over_jterm2_s;
            }
            sum_over_q += g_2pq * (c1 * sum_over_kterm1 + sum_over_kterm2);
            sum_over_q_s += g_2pq * (c1 * sum_over_kterm1_s + sum_over_kterm2_s);
        }
        sum_over_p += f_21p * sum_over_q;
        sum_over_p_s += f_21p * sum_over_q_s;
    }
    let term_m1 = -4. / sqrt!(15_f64) * sum_over_p;
    let term_m1_s = -4. / sqrt!(15_f64) * sum_over_p_s;

    sum_over_p = 0.;
    sum_over_p_s = 0.;
    for (p, (f_22p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[2],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let mut sum_over_q = 0.;
        let mut sum_over_q_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let tmp_q = q_min + q;
            let rek2_22pq = kaula.love_numbers.real(2, p, tmp_q);
            let imk2_22pq = kaula.love_numbers.imaginary(2, p, tmp_q);

            let mut sum_over_kterm1 = 0.;
            let mut sum_over_kterm2 = 0.;
            let mut sum_over_kterm1_s = 0.;
            let mut sum_over_kterm2_s = 0.;

            for (k, (f_31k, f_33k, sum_g_3kj)) in izip!(
                kaula.polynomials.inclination_function_f_3mp[1],
                kaula.polynomials.inclination_function_f_3mp[3],
                kaula.polynomials.eccentricity_function_g_3pq
            )
            .enumerate()
            {
                // Sum over p for each m
                let mut sum_over_jterm1 = 0.;
                let mut sum_over_jterm2 = 0.;
                let mut sum_over_jterm1_s = 0.;
                let mut sum_over_jterm2_s = 0.;

                for (j, g_3kj) in sum_g_3kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_alpha_3 = compute_phase_alpha_3(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );
                    let phase_alpha_4 = compute_phase_alpha_4(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        spin_angle,
                        argument_perihelion,
                        longitude_of_ascending_node,
                        heliocentric_varphi,
                    );

                    let cos_alpha_1 = cos!(phase_alpha_3);
                    let sin_alpha_1 = sin!(phase_alpha_3);
                    let cos_alpha_2 = cos!(phase_alpha_4);
                    let sin_alpha_2 = sin!(phase_alpha_4);

                    sum_over_jterm1 += g_3kj * (sin_alpha_1 * rek2_22pq - cos_alpha_1 * imk2_22pq);
                    sum_over_jterm2 += g_3kj * (sin_alpha_2 * rek2_22pq - cos_alpha_2 * imk2_22pq);
                    sum_over_jterm1_s -= g_3kj * imk2_22pq;
                    sum_over_jterm2_s -= g_3kj * imk2_22pq;
                }
                sum_over_kterm1 += f_33k * sum_over_jterm1;
                sum_over_kterm2 += f_31k * sum_over_jterm2;
                sum_over_kterm1_s += f_33k * sum_over_jterm1_s;
                sum_over_kterm2_s += f_31k * sum_over_jterm2_s;
            }
            sum_over_q += g_2pq * (sum_over_kterm1 + 2. * sum_over_kterm2);
            sum_over_q_s += g_2pq * (sum_over_kterm1_s + 2. * sum_over_kterm2_s);
        }
        sum_over_p += f_22p * sum_over_q;
        sum_over_p_s += f_22p * sum_over_q_s;
    }
    let term_m2 = -5. / (48. * sqrt!(6_f64)) * sum_over_p;
    let term_m2_s = -5. / (48. * sqrt!(6_f64)) * sum_over_p_s;

    let cste_3d = if central_body {
        // Stellar tide: `particle` is the STAR (deformed body), `tidal_host_particle` the PLANET.
        calculate_base_constant(particle, tidal_host_particle, 1.0, semi_major_axis)
    } else {
        // Planetary tide: `particle` is the PLANET (deformed body), `tidal_host_particle` the STAR.
        calculate_base_constant(particle, tidal_host_particle, 1.0, semi_major_axis)
    };

    let orthogonal_force = cste_3d * (term_m1 + term_m2);
    let orthogonal_force_secular = cste_3d * (term_m1_s + term_m2_s);
    (orthogonal_force, orthogonal_force_secular)
}

// radial (e_{r}) component of tidal force
fn calculate_radial_component_of_the_tidal_force_3d(
    tidal_host_particle: &Particle,
    particle: &Particle,
    keplerian_elements: &KeplerianElements,
    central_body: bool,
) -> (f64, f64) {
    let kaula = particle.tides.get_kaula();
    // keplerian elements
    let (
        semi_major_axis,
        _perihelion_distance,
        eccentricity,
        _inclination,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        _orbital_period,
    ) = keplerian_elements.unpack();

    let argument_perihelion = longitude_perihelion - longitude_of_ascending_node;
    // if there is a non-negligible inclination (3D case)
    // Modification has been made but no testing has been done for 3D case (26/2)
    // let eccentricity_function_g_2pq = calculate_eccentricity_function_g_2pq(eccentricity);
    // let inclination_function_f_2mp = calculate_inclination_function_f_2mp(obliquity);
    let cste = if central_body {
        // Stellar tide: `particle` is the STAR (deformed body, radius^5) and
        // `tidal_host_particle` the PLANET (perturber, mass^2, heliocentric distance).
        calculate_base_constant(
            particle,            // deformed body: star
            tidal_host_particle, // perturber: planet
            tidal_host_particle.heliocentric_distance,
            semi_major_axis,
        )
    } else {
        // Planetary tide: `particle` is the PLANET (deformed body), `tidal_host_particle` the STAR.
        calculate_base_constant(
            particle,
            tidal_host_particle,
            particle.heliocentric_distance,
            semi_major_axis,
        )
    };

    let (q_min, q_max) = select_eccentricty_order_q(eccentricity);
    let mut sum_over_p_m0 = 0.;
    let mut sum_over_p_m1 = 0.;
    let mut sum_over_p_m2 = 0.;
    let mut sum_over_p_m0_s = 0.;
    let mut sum_over_p_m1_s = 0.;
    let mut sum_over_p_m2_s = 0.;
    for (p, (f_20p, f_21p, f_22p, sum_g_2pq)) in izip!(
        kaula.polynomials.inclination_function_f_2mp[0],
        kaula.polynomials.inclination_function_f_2mp[1],
        kaula.polynomials.inclination_function_f_2mp[2],
        kaula.polynomials.eccentricity_function_g_2pq
    )
    .enumerate()
    {
        let f_20p_2 = f_20p.powi(2);
        let f_21p_2 = f_21p.powi(2);
        let f_22p_2 = f_22p.powi(2);

        let mut sum_over_q_m0 = 0.;
        let mut sum_over_q_m1 = 0.;
        let mut sum_over_q_m2 = 0.;
        let mut sum_over_q_m0_s = 0.;
        let mut sum_over_q_m1_s = 0.;
        let mut sum_over_q_m2_s = 0.;

        for (q, g_2pq) in sum_g_2pq.iter().take(q_max).skip(q_min).enumerate() {
            let g_2pq_2 = g_2pq.powi(2);

            let tmp_q = q_min + q;
            let rek2_20pq = kaula.love_numbers.real(0, p, tmp_q);
            let imk2_20pq = kaula.love_numbers.imaginary(0, p, tmp_q);
            let rek2_21pq = kaula.love_numbers.real(1, p, tmp_q);
            let imk2_21pq = kaula.love_numbers.imaginary(1, p, tmp_q);
            let rek2_22pq = kaula.love_numbers.real(2, p, tmp_q);
            let imk2_22pq = kaula.love_numbers.imaginary(2, p, tmp_q);

            let mut sum_over_k_m0 = 0.;
            let mut sum_over_k_m1 = 0.;
            let mut sum_over_k_m2 = 0.;

            for (k, (f_20k, f_21k, f_22k, sum_g_2kj)) in izip!(
                kaula.polynomials.inclination_function_f_2mp[0],
                kaula.polynomials.inclination_function_f_2mp[1],
                kaula.polynomials.inclination_function_f_2mp[2],
                kaula.polynomials.eccentricity_function_g_2pq
            )
            .enumerate()
            {
                let mut sum_over_j_m0 = 0.;
                let mut sum_over_j_m1 = 0.;
                let mut sum_over_j_m2 = 0.;
                for (j, g_2kj) in sum_g_2kj.iter().take(q_max).skip(q_min).enumerate() {
                    let phase_alpha = alpha_pqkj(
                        p,
                        q_min + q,
                        k,
                        q_min + j,
                        mean_anomaly,
                        argument_perihelion,
                    );

                    let cos_alpha = cos!(phase_alpha);
                    let sin_alpha = sin!(phase_alpha);

                    sum_over_j_m0 += g_2kj * (cos_alpha * rek2_20pq - sin_alpha * imk2_20pq);
                    sum_over_j_m1 += g_2kj * (cos_alpha * rek2_21pq - sin_alpha * imk2_21pq);
                    sum_over_j_m2 += g_2kj * (cos_alpha * rek2_22pq - sin_alpha * imk2_22pq);
                }
                sum_over_k_m0 += f_20k * sum_over_j_m0;
                sum_over_k_m1 += f_21k * sum_over_j_m1;
                sum_over_k_m2 += f_22k * sum_over_j_m2;
            }
            sum_over_q_m0 += g_2pq * sum_over_k_m0;
            sum_over_q_m1 += g_2pq * sum_over_k_m1;
            sum_over_q_m2 += g_2pq * sum_over_k_m2;
            sum_over_q_m0_s += g_2pq_2 * rek2_20pq;
            sum_over_q_m1_s += g_2pq_2 * rek2_21pq;
            sum_over_q_m2_s += g_2pq_2 * rek2_22pq;
        }
        sum_over_p_m0 += f_20p * sum_over_q_m0;
        sum_over_p_m1 += f_21p * sum_over_q_m1;
        sum_over_p_m2 += f_22p * sum_over_q_m2;
        sum_over_p_m0_s += f_20p_2 * sum_over_q_m0_s;
        sum_over_p_m1_s += f_21p_2 * sum_over_q_m1_s;
        sum_over_p_m2_s += f_22p_2 * sum_over_q_m2_s;
    }

    let radial_force =
        (sum_over_p_m0 + (1. / 3.) * sum_over_p_m1 + (1. / 12.) * sum_over_p_m2) * cste;
    let radial_force_secular =
        (sum_over_p_m0_s + (1. / 3.) * sum_over_p_m1_s + (1. / 12.) * sum_over_p_m2_s) * cste;
    (radial_force, radial_force_secular)
}
