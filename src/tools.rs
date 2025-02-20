use super::constants::{AU, G_SI, HOUR, M_SUN, PI, TWO_PI};
use super::particles::Axes;
use std;
use time;
use time::{OffsetDateTime, format_description};

pub fn calculate_inclination_orbital_equatorial_plane(
    position: Axes,
    velocity: Axes,
    spin: Axes,
) -> f64 {
    // Calculate the spin axis inclination described as
    // the inclination of the orbit plane with respect to the equatorial plane

    let (x, y, z) = position.unpack();
    let (u, v, w) = velocity.unpack();
    let (sx, sy, sz) = spin.unpack();
    let s = sqrt!(sx.powi(2) + sy.powi(2) + sz.powi(2));

    // Calculate the component of the orbital angular momentum
    let hx = y * w - z * v;
    let hy = z * u - x * w;
    let hz = x * v - y * u;
    let h = sqrt!(hx.powi(2) + hy.powi(2) + hz.powi(2));

    let hx_rel = hx / h;
    let hy_rel = hy / h;
    let hz_rel = hz / h;
    let h_rel = sqrt!(hx_rel.powi(2) + hy_rel.powi(2) + hz_rel.powi(2));

    let numerator = hx_rel * sx + hy_rel * sy + hz_rel * sz;
    let denominator = h_rel * s;

    let cos_inclination = (numerator / denominator).clamp(-1., 1.);

    acos!(cos_inclination)
}

pub fn calculate_eccentricity_vector(gm: f64, position: Axes, velocity: Axes) -> Axes {
    let (x, y, z) = position.unpack();
    let (u, v, w) = velocity.unpack();

    // Angular momentum
    let hx = y * w - z * v;
    let hy = z * u - x * w;
    let hz = x * v - y * u;

    // v vectorial h
    let v_vect_h_x = v * hz - w * hy;
    let v_vect_h_y = w * hx - u * hz;
    let v_vect_h_z = u * hy - v * hx;

    // distance
    let r = sqrt!(x * x + y * y + z * z);

    // eccentricity component
    Axes::from(
        (v_vect_h_x / gm) - (x / r),
        (v_vect_h_y / gm) - (y / r),
        (v_vect_h_z / gm) - (z / r),
    )
}

pub fn calculate_keplerian_orbital_elements(
    gm: f64,
    position: Axes,
    velocity: Axes,
) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    // Based on the implementation of Chambers in Mercury
    // Calculates Keplerian orbital elements given relative coordinates and
    // velocities, and GM = G times the sum of the masses.

    let (x, y, z) = position.unpack();
    let (u, v, w) = velocity.unpack();

    // --- Output
    // semi-major axis (in AU)
    // perihelion distance
    let mut eccentricity; // eccentricity
    let mut i = 0.; // inclination
    // longitude of perihelion (NOT argument of perihelion!!)
    // longitude of ascending node
    // mean anomaly (or mean longitude if eccentricity < 1.e-8)
    // orbital_period (in days)

    // --- Local
    let (ex, ey, ez) = calculate_eccentricity_vector(gm, position, velocity).unpack();
    let e_scal_r = ex * x + ey * y + ez * z;
    let hx = y * w - z * v;
    let hy = z * u - x * w;
    let hz = x * v - y * u;
    let h2 = hx.powi(2) + hy.powi(2) + hz.powi(2);
    let v2 = u * u + v * v + w * w;
    let rv = x * u + y * v + z * w;
    let r = sqrt!(x * x + y * y + z * z);
    let h = sqrt!(h2);
    let s = h2 / gm;

    // --- Semi-major axis
    let a = gm * r / (2.0 * gm - r * v2);

    // --- Inclination and longitude of ascending node
    let mut longitude_of_ascending_node;
    let cos_i = (hz / h).clamp(-1., 1.);
    if abs!(cos_i) < 1. {
        i = acos!(cos_i);
        if hy > 0. {
            i = TWO_PI - i;
        }
        //longitude_of_ascending_node = atan!(hx / -hy);
        longitude_of_ascending_node = hx.atan2(-hy);
        if longitude_of_ascending_node < 0. {
            longitude_of_ascending_node += TWO_PI;
        }
        if i == 0. {
            longitude_of_ascending_node = 0.;
        }
    } else {
        if cos_i > 0. {
            i = 0.;
        }
        if cos_i < 0. {
            i = PI;
        }
        longitude_of_ascending_node = 0.;
    }
    //

    // --- Eccentricity and perihelion distance
    let temp = 1. + s * (v2 / gm - 2. / r);
    if temp <= 0. {
        eccentricity = 0.;
    } else {
        eccentricity = sqrt!(temp);
    }
    let q = s / (1. + eccentricity);
    if eccentricity < 3.0e-8 {
        eccentricity = 0.;
    }

    // --- Longitude of perihelion
    //let longitude_perihelion:f64;
    //if eccentricity == 0. {
    //longitude_perihelion= 0.;
    //} else {
    //longitude_perihelion = modulus(ey.atan2(ex) + TWO_PI + TWO_PI, TWO_PI);
    //}

    let argument_perihelion = if eccentricity == 0. {
        0.
    } else {
        let ex_rot =
            ex * cos!(longitude_of_ascending_node) + ey * sin!(longitude_of_ascending_node);
        let ey_rot =
            -ex * sin!(longitude_of_ascending_node) + ey * cos!(longitude_of_ascending_node);
        let ez_rot = ez;
        let ex_rot_2 = ex_rot;
        let ey_rot_2 = ey_rot * cos!(i) + ez_rot * sin!(i);
        //let ez_rot_2 = -ey_rot * sin!(i) + ez_rot * cos!(i);

        // argument_perihelion = atan!(ey_rot_2/ex_rot_2);
        atan2!(ey_rot_2, ex_rot_2)
    };

    let longitude_perihelion = argument_perihelion + longitude_of_ascending_node;

    // --- True anomaly
    let mut true_anomaly;
    let cos_f;
    if eccentricity == 0. {
        let x_rot = x * cos!(longitude_of_ascending_node) + y * sin!(longitude_of_ascending_node);
        cos_f = (x_rot / r).clamp(-1., 1.);
        true_anomaly = acos!(cos_f);
    } else {
        cos_f = (e_scal_r / (eccentricity * r)).clamp(-1., 1.);
        true_anomaly = acos!(cos_f);
        if rv < 0. {
            true_anomaly = TWO_PI - true_anomaly;
        }
    }

    // --- Mean anomaly
    let mut mean_anomaly;
    if eccentricity == 0. {
        mean_anomaly = true_anomaly;
    } else {
        let mut cos_bige = ((1. / eccentricity) * (1. - (r / a))).clamp(-1., 1.);

        // Mean anomaly for ellipse
        if eccentricity < 1. {
            if abs!(cos_bige) > 1. {
                cos_bige = cos_bige.signum();
            }
            let mut bige = acos!(cos_bige);
            if rv < 0. {
                bige = TWO_PI - bige;
            }
            mean_anomaly = bige - eccentricity * sin!(bige);
        } else {
            // Mean anomaly for hyperbola
            if cos_bige < 1. {
                cos_bige = 1.;
            }
            let mut bige = ln!(sqrt!(cos_bige + (cos_bige * cos_bige - 1.)));
            if rv < 0. {
                bige = -bige;
            }
            mean_anomaly = eccentricity * sinh!(bige) - bige;
        }
    }

    if mean_anomaly < 0. {
        mean_anomaly += TWO_PI;
    }
    if mean_anomaly > TWO_PI {
        mean_anomaly = modulus(mean_anomaly, TWO_PI);
    }

    // Given relative coordinates and velocities (of the body 1 respect to body 2),
    // and GM = G times the sum of the masses (body 1 + body 2)
    let orbital_period = (TWO_PI / sqrt!(gm)) * a.powf(3. / 2.); // in days

    // a Semimajor axis
    // q perihelion distance
    // e eccentricity
    // i inclination
    // p longitude of perihelion (! not the argument of perihelion)
    // n longitude of ascending node
    // l mean anomaly
    // orbital period bah...

    (
        a,
        q,
        eccentricity,
        i,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        orbital_period,
    )
}

pub fn calculate_perihelion_distance_and_eccentricity(
    gm: f64,
    position: Axes,
    velocity: Axes,
) -> (f64, f64) {
    //Based on the implementation of Chambers in Mercury
    // Calculates Keplerian orbital elements given relative coordinates and
    // velocities, and GM = G times the sum of the masses.

    let (x, y, z) = position.unpack();
    let (u, v, w) = velocity.unpack();

    // perihelion distance
    // Local
    let hx = y * w - z * v;
    let hy = z * u - x * w;
    let hz = x * v - y * u;
    let h2 = hx.powi(2) + hy.powi(2) + hz.powi(2);
    let v2 = u * u + v * v + w * w;
    let r = sqrt!(x * x + y * y + z * z);
    let s = h2 / gm;

    // Eccentricity and perihelion distance
    let temp = 1. + s * (v2 / gm - 2. / r);

    let eccentricity = if temp <= 0. { 0. } else { sqrt!(temp) };
    let q = s / (1. + eccentricity);

    (q, eccentricity)
}

pub fn calculate_cartesian_coordinates(
    gm: f64,
    q: f64,
    e: f64,
    i0: f64,
    p: f64,
    n0: f64,
    l: f64,
) -> (f64, f64, f64, f64, f64, f64) {
    // Calculates Cartesian coordinates and velocities given Keplerian orbital
    // elements (for elliptical, parabolic or hyperbolic orbits).
    // Based on the implementation of Chambers in Mercury, which is
    //   based on a routine from Levison and Duncan's SWIFT integrator.
    // WARNING: It gives NaN velocities when eccentricity == 1. (also in the original implementation in mercury code)

    // Input
    // gm  = grav const * (central + secondary mass)
    // q  = perihelion distance
    // e  = eccentricity
    // i  = inclination (degrees)
    // p  = longitude of perihelion !!!
    // n  = longitude of the ascending node (degrees)
    // l  = mean anomaly

    // Output
    // Cartesian positions  ( units the same as a )
    // Cartesian positions  ( units the same as a )
    // Cartesian positions  ( units the same as a )
    // Cartesian velocities ( units the same as sqrt(gm/a) )
    // Cartesian velocities ( units the same as sqrt(gm/a) )
    // Cartesian velocities ( units the same as sqrt(gm/a) )

    // Change from longitude of perihelion to argument of perihelion
    let g0 = p - n0;

    // Rotation factors
    let si = sin!(i0);
    let ci = cos!(i0);
    let sg = sin!(g0);
    let cg = cos!(g0);
    let sn = sin!(n0);
    let cn = cos!(n0);
    //let (i, si, ci) = mco_sine(i0);
    //let (g, sg, cg) = mco_sine(g0);
    //let (n, sn, cn) = mco_sine(n0);
    let mut z1 = cg * cn;
    let mut z2 = cg * sn;
    let mut z3 = sg * cn;
    let mut z4 = sg * sn;
    let d11 = z1 - z4 * ci;
    let d12 = z2 + z3 * ci;
    let d13 = sg * si;
    let d21 = -z3 - z2 * ci;
    let d22 = -z4 + z1 * ci;
    let d23 = cg * si;

    // Semi-major axis
    let a = q / (1. - e);

    let se;
    let ce;

    // Ellipse
    if e < 1. {
        let romes = sqrt!(1. - e * e);
        let temp0 = kepler_solution_for_eccentrities_smaller_than_one(e, l);
        se = sin!(temp0);
        ce = cos!(temp0);
        //let (_, se, ce) = mco_sine(temp0);
        z1 = a * (ce - e);
        z2 = a * romes * se;
        let temp = sqrt!(gm / a) / (1. - e * ce);
        z3 = -se * temp;
        z4 = romes * ce * temp;
    } else {
        // Parabola
        if e == 1. {
            let tmp = kepler_solution_for_a_parabola(l);
            ce = tmp.1;
            z1 = q * (1. - ce * ce);
            z2 = 2. * q * ce;
            z4 = sqrt!(2. * gm / q) / (1. + ce * ce);
            z3 = -ce * z4;
        } else {
            // Hyperbola
            let romes = sqrt!(e * e - 1.);
            let mut temp = kepler_solution_for_a_hyperbola(e, l);
            se = sinh!(temp);
            ce = cosh!(temp);
            z1 = a * (ce - e);
            z2 = -a * romes * se;
            temp = sqrt!(gm / abs!(a)) / (e * ce - 1.);
            z3 = -se * temp;
            z4 = romes * ce * temp;
        }
    }

    let x = d11 * z1 + d21 * z2;
    let y = d12 * z1 + d22 * z2;
    let z = d13 * z1 + d23 * z2;
    let u = d11 * z3 + d21 * z4;
    let v = d12 * z3 + d22 * z4;
    let w = d13 * z3 + d23 * z4;

    (x, y, z, u, v, w)
}

fn kepler_solution_for_eccentrities_smaller_than_one(e: f64, oldl: f64) -> f64 {
    //Based on the implementation of Chambers in Mercury (mco_kep)
    // Solves Kepler's equation for eccentricities less than one.
    // Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
    //
    // e = eccentricity\n
    // l = mean anomaly      (radians)\n
    // u = eccentric anomaly (   "   )\n

    //// Local
    //real(double_precision) :: l,pi,twopi,piby2,u1,u2,ome,sign
    //real(double_precision) :: x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
    //real(double_precision) :: p,q,p2,ss,cc
    //logical flag,big,bigg

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
    let q;
    let sn;

    let mut l;
    let piby2 = 0.5 * PI;

    // Reduce mean anomaly to lie in the range 0 < l < pi
    if oldl >= 0. {
        l = modulus(oldl, TWO_PI);
    } else {
        l = (modulus(oldl, TWO_PI)) + TWO_PI;
    }
    let mut sign = 1.;
    if l > PI {
        l = TWO_PI - l;
        sign = -1.;
    }

    let ome = 1. - e;

    if (l >= 0.45) || (e < 0.55) {
        // Regions A,B or C in Nijenhuis
        // -----------------------------

        // Rough starting value for eccentric anomaly
        if l < ome {
            u1 = ome;
        } else if l > (PI - 1. - e) {
            u1 = (l + e * PI) / (1. + e);
        } else {
            u1 = l + e;
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
        f2 = e * sn;
        f0 = u1 - f2 - l;
        f1 = 1. - e * dsn;
        u2 = u1 - f0 / (f1 - 0.5 * f0 * f2 / f1);
    } else {
        // Region D in Nijenhuis
        // ---------------------

        // Rough starting value for eccentric anomaly
        z1 = 4. * e + 0.5;
        p = ome / z1;
        q = 0.5 * l / z1;
        p2 = p * p;
        z2 = ((sqrt!(p2 * p + q * q) + q).log(std::f64::consts::E) / 1.5).exp();
        u1 = 2. * q / (z2 + p + p2 / z2);

        // Improved value using Newton's method
        z2 = u1 * u1;
        z3 = z2 * z2;
        u2 = u1 - 0.075 * u1 * z3 / (ome + z1 * z2 + 0.375 * z3);
        u2 = l + e * u2 * (3. - 4. * u2 * u2);
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

    f0 = l - u2 * ome - e * z1;
    f1 = ome + e * z2;
    f2 = 0.5 * e * (u2 - z1);
    let f3 = e / 6. * (1. - z2);
    z1 = f0 / f1;
    z2 = f0 / (f2 * z1 + f1);
    sign * (u2 + f0 / ((f3 * z1 + f2) * z2 + f1))
}

fn kepler_solution_for_a_parabola(mut q: f64) -> (f64, f64) {
    // Based on the implementation of Duncan in Mercury (orbel_zget)
    // Solves the equivalent of Kepler's eqn. for a parabola
    // given the parabola mean anomaly (q) (Fitz. notation.)
    //
    // ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
    //
    // @remarks For a parabola we can solve analytically.

    let inverted = if q < 0. {
        q = abs!(q);
        true
    } else {
        false
    };

    let eccentric_anomaly = if q < 1.0e-3 {
        q * (1. - (q * q / 3.) * (1. - q * q))
    } else {
        let tmp = (0.5 * (3. * q + sqrt!(9. * q.powi(2) + 4.))).powf(1. / 3.);
        tmp - 1. / tmp
    };

    if inverted {
        (-q, -eccentric_anomaly)
    } else {
        (q, eccentric_anomaly)
    }
}

fn kepler_solution_for_a_hyperbola(eccentricity: f64, hyperbola_mean_anomaly: f64) -> f64 {
    // Based on the implementation of Duncan in Mercury (orbel_fhybrid)
    // Solves Kepler's eqn. for hyperbola using hybrid approach.

    // orbel_fhybrid (eccentric anomaly)
    if abs!(hyperbola_mean_anomaly) < 0.636 * eccentricity - 0.6 {
        let tmp = kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(
            eccentricity,
            hyperbola_mean_anomaly,
        );
        tmp.1
    } else {
        kepler_solution_for_a_hyperbola_hybrid_approach(eccentricity, hyperbola_mean_anomaly)
    }
}

fn kepler_solution_for_a_hyperbola_hybrid_approach_for_low_n(e: f64, capn0: f64) -> (f64, f64) {
    //Based on the implementation of Duncan in Mercury (orbel_flon)
    // Solves Kepler's eqn. for hyperbola using hybrid approach. \n\n
    // ALGORITHM: Uses power series for N in terms of F and Newton,s method
    // REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)

    // Inputs:
    // e = eccentricty
    // capn = hyperbola mean anomaly

    // Output
    let mut orbel_flon;

    // copy of the input capn0 that is not modified
    let mut capn = capn0;
    let mut iflag: i32;

    let mut x;
    let mut x2;
    let mut f;
    let mut fp;
    let mut dx;
    let diff;

    let imax: i32 = 10;
    let tiny = 4.0e-15; // A small number
    let a3 = 1_037_836_800.;
    let a5 = 51_891_840.;
    let a7 = 1_235_520.;
    let a9 = 17_160.;
    let a11 = 156.;

    let b3 = 3. * a3;
    let b5 = 5. * a5;
    let b7 = 7. * a7;
    let b9 = 9. * a9;
    let b11 = 11. * a11;

    // Function to solve "Kepler's eqn" for F (here called
    // x) for given e and CAPN. Only good for smallish CAPN

    iflag = 0;
    if capn < 0. {
        iflag = 1;
        capn = -capn;
    }

    let a_tmp = 6_227_020_800.;
    let a1 = a_tmp * (1. - 1. / e);
    let a0 = -a_tmp * capn / e;
    let b1 = a1;

    //  Set iflag nonzero if capn < 0., in which case solve for -capn
    //  and change the sign of the final answer for F.
    //  Begin with a reasonable guess based on solving the cubic for small F

    let a = 6. * (e - 1.) / e;
    let b = -6. * capn / e;
    let sq = sqrt!(0.25 * b * b + a * a * a / 27.);
    let biga = (-0.5 * b + sq).powf(1. / 3.);
    let bigb = -(0.5 * b + sq).powf(1. / 3.);
    x = biga + bigb;
    //  write(6,*) 'cubic = ',x**3 +a*x +b
    orbel_flon = x;
    // If capn is tiny (or zero) no need to go further than cubic even for
    // e =1.
    if capn >= tiny {
        let mut converge: bool = false;
        for _ in 1..imax {
            x2 = x * x;
            f = a0 + x * (a1 + x2 * (a3 + x2 * (a5 + x2 * (a7 + x2 * (a9 + x2 * (a11 + x2))))));
            fp = b1 + x2 * (b3 + x2 * (b5 + x2 * (b7 + x2 * (b9 + x2 * (b11 + 13. * x2)))));
            dx = -f / fp;
            orbel_flon = x + dx;
            //   If we have converged here there's no point in going on
            if abs!(dx) < tiny {
                converge = true;
                break;
            }
            x = orbel_flon;
        }

        if !converge {
            // Abnormal return here - we've gone thru the loop
            // IMAX times without convergence
            if iflag == 1 {
                orbel_flon = -orbel_flon;
                capn = -capn;
            }
            println!(
                "[WARNING {} UTC] FLON : RETURNING WITHOUT COMPLETE CONVERGENCE",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            diff = e * sinh!(orbel_flon) - orbel_flon - capn;
            println!("N, F, ecc * F.sinh() - F - N : ");
            println!("{capn} {orbel_flon} {diff}");
            return (orbel_flon, capn);
        }
    }

    //  Normal return here, but check if capn was originally negative
    if iflag == 1 {
        orbel_flon = -orbel_flon;
        capn = -capn;
    }

    (orbel_flon, capn)
}

fn kepler_solution_for_a_hyperbola_hybrid_approach(e: f64, capn: f64) -> f64 {
    //Based on the implementation of Duncan in Mercury (orbel_fget)
    // Solves Kepler's eqn. for hyperbola using hybrid approach.
    // ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
    //              Cel. Mech. ".  Quartic convergence from Danby's book.

    // Inputs:
    // e = eccentricty
    // capn = hyperbola mean anomaly

    // Output
    let mut orbel_fget; // eccentric anomaly

    //...  Internals:
    let imax: i32 = 10;
    let tiny = 4.0e-15; // A small number

    let tmp;
    let mut x;
    let mut shx;
    let mut chx;
    let mut esh;
    let mut ech;
    let mut f;
    let mut fp;
    let mut fpp;
    let mut fppp;
    let mut dx;

    //----
    // Function to solve "Kepler's eqn" for F (here called
    // x) for given e and CAPN.

    //  begin with a guess proposed by Danby
    if capn < 0. {
        tmp = -2. * capn / e + 1.8;
        x = -tmp.log(std::f64::consts::E);
    } else {
        tmp = 2. * capn / e + 1.8;
        x = tmp.log(std::f64::consts::E);
    }

    orbel_fget = x;

    for _ in 1..imax {
        shx = sinh!(x);
        chx = cosh!(x);
        esh = e * shx;
        ech = e * chx;
        f = esh - x - capn;
        fp = ech - 1.;
        fpp = esh;
        fppp = ech;
        dx = -f / fp;
        dx = -f / (fp + dx * fpp / 2.);
        dx = -f / (fp + dx * fpp / 2. + dx * dx * fppp / 6.);
        orbel_fget = x + dx;
        //   If we have converged here there's no point in going on
        if abs!(dx) <= tiny {
            return orbel_fget;
        }
        x = orbel_fget;
    }

    println!(
        "[WARNING {} UTC] FGET : RETURNING WITHOUT COMPLETE CONVERGENCE",
        OffsetDateTime::now_utc()
            .format(
                &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                    .unwrap()
            )
            .unwrap()
    );
    orbel_fget
}

#[allow(dead_code)]
fn mco_sine(x0: f64) -> (f64, f64, f64) {
    //Based on the implementation of Chambers in Mercury
    // Calculates sin and cos of an angle X (in radians).
    // Chambers says: "x is modified by the routine. Without this part, outputs of mercury are different, don't ask me why."

    // TODO why results of this routine are different from simple calls of intrinsec cos() and sin()
    let x = if x0 > 0. {
        modulus(x0, TWO_PI)
    } else {
        modulus(x0, TWO_PI + TWO_PI)
    };

    let cx = x.cos();

    let sx = if x > PI {
        -sqrt!(1.0 - cx.powi(2))
    } else {
        sqrt!(1.0 - cx.powi(2))
    };

    (x, sx, cx)
}

fn modulus(a: f64, b: f64) -> f64 {
    //println!("Modulus: {}", -21.0f64 % 4.0f64);         // -1 because -21 divided by 4 gives -5 with a remainder of -1.
    //println!("Modulus: {}", 21.0f64 % 4.0f64);          //  1
    //println!("Modulus: {}", modulus(-21.0f64, 4.0f64)); //  3 because -21 + 4 x 6 is 3.
    //println!("Modulus: {}", modulus(21.0f64, 4.0f64));  //  1
    a - (a / b).floor() * b
}

pub fn linear_interpolation(target_x: f64, x: &[f64], y: &[f64]) -> (f64, usize) {
    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to target_x.
    let (left, right) = find_indices_around_target_value(x, target_x);

    let target_y = if left == right {
        // Target value out of range, use limit values
        y[left]
    } else {
        //// Interpolate
        // Linear
        //target_y = (y[left] * (x[right] - target_x) + y[right] * (target_x - x[left])) / (x[right] - x[left])
        // Linear (alternative)
        let x_left = x[left];
        let target_percent = (target_x - x_left) / (x[right] - x_left); // Transform target to percent as in transforming x[left]..x[right] to 0..1
        y[left] * (1. - target_percent) + y[right] * target_percent
    };

    (target_y, left)
}

pub fn cosine_interpolation(target_x: f64, x: &[f64], y: &[f64]) -> (f64, usize) {
    let target_y;

    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to target_x.
    let (left, right) = find_indices_around_target_value(x, target_x);

    if left == right {
        // Target value out of range, use limit values
        target_y = y[left];
    } else {
        //// Interpolate
        // Cosine interpolate (http://paulbourke.net/miscellaneous/interpolation/)
        // - Smooth around the real data points (contrary to the linear interpolation)
        let x_left = x[left];
        let mut target_percent = (target_x - x_left) / (x[right] - x_left); // Transform target to percent as in transforming x[left]..x[right] to 0..1
        target_percent = (1. - cos!(target_percent * PI)) / 2.; // Transform target percent so that it gets smoothed when close to 0 or 1 (i.e., closer to x[left] or x[right])
        target_y = y[left] * (1. - target_percent) + y[right] * target_percent;
    }

    (target_y, left)
}

fn find_indices_around_target_value(data: &[f64], target_value: f64) -> (usize, usize) {
    // Find the nearest interval [ x(LEFT), x(RIGHT) ] to XVAL.
    let ndata = data.len();
    let last_idx = ndata.wrapping_sub(1);
    let (left_idx, right_idx) = match data.iter().position(|&r| r > target_value) {
        None => {
            if data[last_idx] > target_value {
                (0, 0)
            } else {
                (last_idx, last_idx)
            }
        }
        Some(i) => {
            if i == 0 {
                (0, 0)
            } else if i == ndata {
                (last_idx, last_idx)
            } else {
                (i - 1, i)
            }
        }
    };
    (left_idx, right_idx)
}

pub fn calculate_pseudo_synchronization_period(
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
    let pseudo_rot = alpha * sqrt!(G_SI * M_SUN * (star_mass + planet_mass));
    let angular_frequency = pseudo_rot * (semi_major_axis * AU).powf(-3. / 2.) * HOUR * 24.; // days^-1
    // days
    TWO_PI / (angular_frequency)
}

pub fn calculate_spin(
    angular_frequency: f64,
    inclination: f64,
    obliquity: f64,
    longitude_ascending_node: f64,
) -> Axes {
    // Spin taking into consideration the inclination:
    Axes::from(
        angular_frequency * sin!(obliquity + inclination) * sin!(longitude_ascending_node),
        -angular_frequency * sin!(obliquity + inclination) * cos!(longitude_ascending_node),
        angular_frequency * cos!(obliquity + inclination),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{DEG2RAD, G, M_EARTH};
    #[test]
    fn calculate_keplerian_orbital_elements() {
        //---- Star (central body)
        let star_mass = 0.08; // Solar masses
        let planet_mass = 1.0 * M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)

        ////////// Specify initial position and velocity for a stable orbit
        ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
        let a = 0.018; // semi-major axis (in AU)
        let e = 0.1; // eccentricity
        let i = 5. * DEG2RAD; // inclination (degrees)
        let mut p = 0.; // argument of pericentre (degrees)
        let n = 0. * DEG2RAD; // longitude of the ascending node (degrees)
        let l = 0. * DEG2RAD; // mean anomaly (degrees)
        p = (p + n) * DEG2RAD; // Convert to longitude of perihelion !!
        let q = a * (1.0 - e); // perihelion distance
        let gm = G * (planet_mass + star_mass);
        let (x, y, z, vx, vy, vz) = calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

        assert_eq!(x, 0.0162);
        assert_eq!(y, 0.000000000000006571920583533768);
        assert_eq!(z, 0.0000000000000005749685486518799);
        assert_eq!(vx, -0.000000000000014818017716591765);
        assert_eq!(vy, 0.03987438104619194);
        assert_eq!(vz, 0.003488556306654768);
    }
}
