use posidonius::Axes;
use posidonius::constants::{AU, G_SI, HOUR, M_SUN, PI, TWO_PI};
use posidonius::tools::find_indices_around_target_value;

use std;
use time;
use time::{OffsetDateTime, format_description};

pub fn calculate_spin(
    angular_frequency: f64,
    inclination: f64,
    obliquity: f64,
    longitude_ascending_node: f64,
) -> Axes {
    // Spin taking into consideration the inclination:
    Axes::from(
        angular_frequency * f64::sin(obliquity + inclination) * f64::sin(longitude_ascending_node),
        -angular_frequency * f64::sin(obliquity + inclination) * f64::cos(longitude_ascending_node),
        angular_frequency * f64::cos(obliquity + inclination),
    )
}

#[allow(dead_code)]
fn calculate_pseudo_synchronization_period(
    semi_major_axis: f64,
    eccentricity: f64,
    star_mass: f64,
    planet_mass: f64,
) -> f64 {
    let alpha = (1.
        + 15. / 2. * eccentricity.powi(2)
        + 45. / 8. * eccentricity.powi(4)
        + 5. / 16. * eccentricity.powi(6))
        * 1.
        / (1. + 3. * eccentricity.powi(2) + 3. / 8. * eccentricity.powi(4))
        * 1.
        / (1. - eccentricity.powi(2)).powf(1.5);
    let pseudo_rot = alpha * f64::sqrt(G_SI * M_SUN * (star_mass + planet_mass));
    let angular_frequency = pseudo_rot * (semi_major_axis * AU).powf(-3. / 2.) * HOUR * 24.; // days^-1
    // days
    TWO_PI / (angular_frequency)
}

#[allow(dead_code)]
fn cosine_interpolation(target_x: f64, x: &[f64], y: &[f64]) -> (f64, usize) {
    let target_y;

    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to target_x.
    let (left, right) = find_indices_around_target_value(x, target_x);

    if left == right {
        // Target value out of range, use limit values
        target_y = y[left];
    } else {
        //// Interpolate
        // Cosine interpolate (http://paulbourke.net/miscellaneous/interpolation/)
        // Smooth around the real data points (contrary to the linear interpolation)
        let x_left = x[left];
        let mut target_percent = (target_x - x_left) / (x[right] - x_left); // Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_percent = (1. - f64::cos(target_percent * PI)) / 2.; // Transform target percent so that it gets smoothed when close to 0 or 1 (i.e., closer to x[left] or x[right])
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent;
    }

    (target_y, left)
}

// TODO fix this function so it doesn't return NaN
// Calculates Cartesian coordinates and velocities given Keplerian
// orbital elements (for elliptical, parabolic or hyperbolic orbits).
// Based on the implementation of Chambers in Mercury, which is
// based on a routine from Levison and Duncan's SWIFT integrator.
// WARNING: It gives NaN velocities when eccentricity == 1. (also in the original implementation in mercury code)
pub fn calculate_cartesian_coordinates(
    gm: f64,
    perihelion_distance: f64,
    eccentricity: f64,
    inclination: f64,
    longitude_perihelion: f64,
    longitude_of_ascending_node: f64,
    mean_anomaly: f64,
) -> (Axes, Axes) {
    // Input
    // gm  = grav const * (central + secondary mass)

    // Change from longitude of perihelion to argument of perihelion
    let argument_of_perihelion = longitude_perihelion - longitude_of_ascending_node;

    // Rotation factors
    let z1 = f64::cos(argument_of_perihelion) * f64::cos(longitude_of_ascending_node);
    let z2 = f64::cos(argument_of_perihelion) * f64::sin(longitude_of_ascending_node);
    let z3 = f64::sin(argument_of_perihelion) * f64::cos(longitude_of_ascending_node);
    let z4 = f64::sin(argument_of_perihelion) * f64::sin(longitude_of_ascending_node);
    let d11 = z1 - z4 * f64::cos(inclination);
    let d12 = z2 + z3 * f64::cos(inclination);
    let d13 = f64::sin(argument_of_perihelion) * f64::sin(inclination);
    let d21 = -z3 - z2 * f64::cos(inclination);
    let d22 = -z4 + z1 * f64::cos(inclination);
    let d23 = f64::cos(argument_of_perihelion) * f64::sin(inclination);

    let semi_major_axis = perihelion_distance / (1. - eccentricity);

    let (z1, z2, z3, z4) = if eccentricity < 1. {
        // Ellipse
        let romes = f64::sqrt(1. - eccentricity.powi(2));
        let temp0 = kepler_solution_for_eccentrities_smaller_than_one(eccentricity, mean_anomaly);
        let z1 = semi_major_axis * (f64::cos(temp0) - eccentricity);
        let z2 = semi_major_axis * romes * f64::sin(temp0);
        let temp = f64::sqrt(gm / semi_major_axis) / (1. - eccentricity * f64::cos(temp0));
        let z3 = -f64::sin(temp0) * temp;
        let z4 = romes * f64::cos(temp0) * temp;

        (z1, z2, z3, z4)
    } else if eccentricity == 1. {
        // Parabola
        let eccentric_anomaly = kepler_solution_for_a_parabola(mean_anomaly);
        let z1 = perihelion_distance * (1. - eccentric_anomaly.powi(2));
        let z2 = 2. * perihelion_distance * eccentric_anomaly;
        let z4 = f64::sqrt(2. * gm / perihelion_distance) / (1. + eccentric_anomaly.powi(2));
        let z3 = -eccentric_anomaly * z4;

        (z1, z2, z3, z4)
    } else {
        // Hyperbola
        let romes = f64::sqrt(eccentricity.powi(2) - 1.);
        let temp = kepler_solution_for_a_hyperbola(eccentricity, mean_anomaly);
        let z1 = semi_major_axis * (f64::cosh(temp) - eccentricity);
        let z2 = -semi_major_axis * romes * f64::sinh(temp);
        let temp =
            f64::sqrt(gm / f64::abs(semi_major_axis)) / (eccentricity * f64::cosh(temp) - 1.);
        let z3 = -f64::sinh(temp) * temp;
        let z4 = romes * f64::cosh(temp) * temp;

        (z1, z2, z3, z4)
    };

    // Cartesian positions  ( units the same as semi_major_axis )
    let positions = Axes::from(
        d11 * z1 + d21 * z2,
        d12 * z1 + d22 * z2,
        d13 * z1 + d23 * z2,
    );

    // Cartesian velocities ( units the same as sqrt(gm/a) )
    let velocities = Axes::from(
        d11 * z3 + d21 * z4,
        d12 * z3 + d22 * z4,
        d13 * z3 + d23 * z4,
    );

    (positions, velocities)
}

// TODO reference
// Based on the implementation of Chambers in Mercury (mco_kep)
// Solves Kepler's equation for eccentricities less than one.
// Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
fn kepler_solution_for_eccentrities_smaller_than_one(
    eccentricity: f64,
    old_mean_anomaly: f64,
) -> f64 {
    // u = eccentric anomaly (   "   )

    let mut x;
    let mut x2;
    let mut dsn;
    let u1;
    let mut u2;
    let mut f0;
    let mut f1;
    let mut f2;

    let mut z1;
    let mut z2;
    let mut z3;
    let p;
    let p2;
    let sn;

    let mut mean_anomaly;
    let piby2 = 0.5 * PI;

    // Reduce mean anomaly to lie in the range 0 < mean_anomaly < pi
    if old_mean_anomaly >= 0. {
        mean_anomaly = old_mean_anomaly.rem_euclid(TWO_PI);
    } else {
        mean_anomaly = old_mean_anomaly.rem_euclid(TWO_PI) + TWO_PI;
    }
    let mut sign = 1.;
    if mean_anomaly > PI {
        mean_anomaly = TWO_PI - mean_anomaly;
        sign = -1.;
    }

    let ome = 1. - eccentricity;

    if (mean_anomaly >= 0.45) || (eccentricity < 0.55) {
        // Regions A,B or C in Nijenhuis

        // Rough starting value for eccentric anomaly
        if mean_anomaly < ome {
            u1 = ome;
        } else if mean_anomaly > (PI - 1. - eccentricity) {
            u1 = (mean_anomaly + eccentricity * PI) / (1. + eccentricity);
        } else {
            u1 = mean_anomaly + eccentricity;
        }

        // Improved value using Halley's method
        let flag = u1 > piby2;
        if flag {
            x = PI - u1;
        } else {
            x = u1;
        }
        x2 = x * x;
        sn = x * (1. + x2 * (-0.16605 + x2 * 0.00761));
        dsn = 1. + x2 * (-0.49815 + x2 * 0.03805);
        if flag {
            dsn = -dsn;
        }
        f2 = eccentricity * sn;
        f0 = u1 - f2 - mean_anomaly;
        f1 = 1. - eccentricity * dsn;
        u2 = u1 - f0 / (f1 - 0.5 * f0 * f2 / f1);
    } else {
        // Region D in Nijenhuis

        // Rough starting value for eccentric anomaly
        z1 = 4. * eccentricity + 0.5;
        p = ome / z1;
        let perihelion_distance = 0.5 * mean_anomaly / z1;
        p2 = p * p;
        z2 = ((f64::sqrt(p2 * p + perihelion_distance.powi(2)) + perihelion_distance)
            .log(std::f64::consts::E)
            / 1.5)
            .exp();
        u1 = 2. * perihelion_distance / (z2 + p + p2 / z2);

        // Improved value using Newton's method
        z2 = u1 * u1;
        z3 = z2 * z2;
        u2 = u1 - 0.075 * u1 * z3 / (ome + z1 * z2 + 0.375 * z3);
        u2 = mean_anomaly + eccentricity * u2 * (3. - 4. * u2 * u2);
    }

    // Accurate value using 3rd-order version of Newton's method
    // N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy

    // First get accurate values for u2 - sin(u2) and 1 - cos(u2)
    let bigg = u2 > piby2;
    if bigg {
        z3 = PI - u2;
    } else {
        z3 = u2;
    }

    let big = z3 > (0.5 * piby2);
    if big {
        x = piby2 - z3;
    } else {
        x = z3;
    }

    x2 = x * x;

    let ss = x * x2 / 6.
        * (1.
            - x2 / 20.
                * (1.
                    - x2 / 42.
                        * (1.
                            - x2 / 72.
                                * (1.
                                    - x2 / 110.
                                        * (1.
                                            - x2 / 156. * (1. - x2 / 210. * (1. - x2 / 272.)))))));
    let cc = x2 / 2.
        * (1.
            - x2 / 12.
                * (1.
                    - x2 / 30.
                        * (1.
                            - x2 / 56.
                                * (1.
                                    - x2 / 90.
                                        * (1.
                                            - x2 / 132.
                                                * (1.
                                                    - x2 / 182.
                                                        * (1. - x2 / 240. * (1. - x2 / 306.))))))));

    if big {
        z1 = cc + z3 - 1.;
        z2 = ss + z3 + 1. - piby2;
    } else {
        z1 = ss;
        z2 = cc;
    }

    if bigg {
        z1 = 2. * u2 + z1 - PI;
        z2 = 2. - z2;
    }

    f0 = mean_anomaly - u2 * ome - eccentricity * z1;
    f1 = ome + eccentricity * z2;
    f2 = 0.5 * eccentricity * (u2 - z1);
    let f3 = eccentricity / 6. * (1. - z2);
    z1 = f0 / f1;
    z2 = f0 / (f2 * z1 + f1);
    sign * (u2 + f0 / ((f3 * z1 + f2) * z2 + f1))
}

// TODO reference
// For a parabola we can solve analytically.
// Based on the implementation of Duncan in Mercury (orbel_zget)
// Solves the equivalent of Kepler's equation for a parabola
// given the parabola mean anomaly (Fitz. notation.)
// p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
fn kepler_solution_for_a_parabola(mut mean_anomaly: f64) -> f64 {
    let sign = mean_anomaly.signum();
    mean_anomaly = f64::abs(mean_anomaly);

    let eccentric_anomaly = if mean_anomaly < 1.0e-3 {
        mean_anomaly * (1. - (mean_anomaly.powi(2) / 3.) * (1. - mean_anomaly.powi(2)))
    } else {
        let tmp =
            (0.5 * (3. * mean_anomaly + f64::sqrt(9. * mean_anomaly.powi(2) + 4.))).powf(1. / 3.);
        tmp - 1. / tmp
    };

    eccentric_anomaly.copysign(sign)
}

// Based on the implementation of Duncan in Mercury (orbel_fhybrid)
// Solves Kepler's equation for hyperbola using hybrid approach.
fn kepler_solution_for_a_hyperbola(eccentricity: f64, hyperbola_mean_anomaly: f64) -> f64 {
    if f64::abs(hyperbola_mean_anomaly) < 0.636 * eccentricity - 0.6 {
        let tmp = kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(
            eccentricity,
            hyperbola_mean_anomaly,
        );
        tmp.1
    } else {
        kepler_solution_for_a_hyperbola_hybrid_approach(eccentricity, hyperbola_mean_anomaly)
    }
}

//Based on the implementation of Duncan in Mercury (orbel_flon)
// Solves Kepler's equation for hyperbola using hybrid approach.
// Uses power series for N in terms of F and Newton's method
// REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
fn kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(
    eccentricity: f64,
    mut hyperbola_mean_anomaly: f64,
) -> (f64, f64) {
    // If `hyperbola_mean_anomaly` is negative, solve for the positive value and reapply the negative sign to the output.
    let sign = hyperbola_mean_anomaly.signum();
    hyperbola_mean_anomaly = f64::abs(hyperbola_mean_anomaly);

    let a_tmp = 6_227_020_800.;
    let a0 = -a_tmp * hyperbola_mean_anomaly / eccentricity;
    let a1 = a_tmp * (1. - 1. / eccentricity);
    let a3 = 1_037_836_800.;
    let a5 = 51_891_840.;
    let a7 = 1_235_520.;
    let a9 = 17_160.;
    let a11 = 156.;

    let b1 = a1;
    let b3 = 3. * a3;
    let b5 = 5. * a5;
    let b7 = 7. * a7;
    let b9 = 9. * a9;
    let b11 = 11. * a11;

    //  Begin with a reasonable guess based on solving the cubic for small F
    let a = 6. * (eccentricity - 1.) / eccentricity;
    let b = -6. * hyperbola_mean_anomaly / eccentricity;
    let sq = f64::sqrt(0.25 * b * b + a * a * a / 27.);
    let biga = (-0.5 * b + sq).powf(1. / 3.);
    let bigb = -(0.5 * b + sq).powf(1. / 3.);
    let mut orbel_flon = biga + bigb;
    let tiny = 4.0e-15; // A small number
    // If hyperbola_mean_anomaly is tiny (or zero) no need to go further than cubic even for eccentricity = 1.
    if hyperbola_mean_anomaly >= tiny {
        let mut converge: bool = false;
        for _ in 1..10 {
            let x2 = orbel_flon.powi(2);
            let f = a0
                + orbel_flon
                    * (a1 + x2 * (a3 + x2 * (a5 + x2 * (a7 + x2 * (a9 + x2 * (a11 + x2))))));
            let fp = b1 + x2 * (b3 + x2 * (b5 + x2 * (b7 + x2 * (b9 + x2 * (b11 + 13. * x2)))));
            let dx = -f / fp;
            orbel_flon += dx;
            //   If we have converged here there's no point in going on
            if f64::abs(dx) < tiny {
                converge = true;
                break;
            }
        }

        if !converge {
            // Abnormal return
            // No convergence after 10 iterations
            println!(
                "[WARNING {} UTC] FLON : RETURNING WITHOUT COMPLETE CONVERGENCE",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            println!(
                "hyperbola_mean_anomaly: {} orbel_flon: {}  sinh(orbel_flon) - orbel_flon - hyperbola_mean_anomaly: {}",
                hyperbola_mean_anomaly.copysign(sign),
                orbel_flon.copysign(sign),
                eccentricity * f64::sinh(orbel_flon.copysign(sign))
                    - orbel_flon.copysign(sign)
                    - hyperbola_mean_anomaly.copysign(sign)
            );
        }
    }

    (
        orbel_flon.copysign(sign),
        hyperbola_mean_anomaly.copysign(sign),
    )
}

// TODO reference
// Based on the implementation of Duncan in Mercury (orbel_fget)
// Solves Kepler's equation for hyperbola using hybrid approach.
// ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
// Cel. Mech. ".  Quartic convergence from Danby's book.
fn kepler_solution_for_a_hyperbola_hybrid_approach(
    eccentricity: f64,
    hyperbola_mean_anomaly: f64,
) -> f64 {
    //  begin with a guess proposed by Danby
    let sign = hyperbola_mean_anomaly.signum();
    let tmp = (2. * hyperbola_mean_anomaly / eccentricity + 1.8).copysign(sign);
    let mut eccentric_anomaly = (f64::ln(tmp)).copysign(sign);

    let mut converged = false;
    for _ in 1..10 {
        let esh = eccentricity * f64::sinh(eccentric_anomaly);
        let ech = eccentricity * f64::cosh(eccentric_anomaly);
        let f = esh - eccentric_anomaly - hyperbola_mean_anomaly;
        let fp = ech - 1.;
        let dx = -f / fp;
        let dx = -f / (fp + dx * esh / 2.);
        let dx = -f / (fp + dx * esh / 2. + dx * dx * ech / 6.);
        eccentric_anomaly += dx;
        // Break as soon as convergence is reached.
        if f64::abs(dx) <= 4.0e-15 {
            // A small number
            converged = true;
            break;
        }
    }

    if !converged {
        println!(
            "[WARNING {} UTC] FGET : RETURNING WITHOUT COMPLETE CONVERGENCE",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
    }

    eccentric_anomaly
}

#[test]
fn calculate_keplerian_orbital_elements() {
    use posidonius::constants::{G, M_EARTH};
    //- Star (central body)
    let star_mass = 0.08; // Solar masses
    let planet_mass = 1.0 * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let semi_major_axis = 0.018; // semi-major axis (in AU)
    let eccentricity = 0.1; // eccentricity
    let inclination = 5_f64.to_radians();
    let mut p = 0.; // argument of pericentre (degrees)
    let n = 0_f64.to_radians(); // longitude of the ascending node (degrees)
    let mean_anomaly = 0_f64.to_radians(); // mean anomaly (degrees)
    p = (p + n).to_radians(); // Convert to longitude of perihelion !!
    let perihelion_distance = semi_major_axis * (1.0 - eccentricity); // perihelion distance
    let gm = G * (planet_mass + star_mass);
    let (position, velocity) = calculate_cartesian_coordinates(
        gm,
        perihelion_distance,
        eccentricity,
        inclination,
        p,
        n,
        mean_anomaly,
    );
    assert_eq!(position.x, 0.0162);
    assert_eq!(position.y, 0.000000000000006571920583533768);
    assert_eq!(position.z, 0.0000000000000005749685486518799);
    assert_eq!(velocity.x, -0.000000000000014818017716591765);
    assert_eq!(velocity.y, 0.03987438104619194);
    assert_eq!(velocity.z, 0.003488556306654768);
}
