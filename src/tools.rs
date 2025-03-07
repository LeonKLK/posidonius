use crate::constants::{PI, TWO_PI};
use crate::particles::Axes;

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

fn calculate_eccentricity_vector(gm: f64, position: Axes, velocity: Axes) -> Axes {
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

pub struct KeplerianElements {
    pub semi_major_axis: f64, // (AU)
    pub perihelion_distance: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub longitude_perihelion: f64, // NOT the argument of perihelion
    pub longitude_of_ascending_node: f64,
    pub mean_anomaly: f64,   // mean longitude if eccentricity < 1.e-8
    pub orbital_period: f64, // (days)
}

impl KeplerianElements {
    pub fn unpack(&self) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
        (
            self.semi_major_axis,
            self.perihelion_distance,
            self.eccentricity,
            self.inclination,
            self.longitude_perihelion,
            self.longitude_of_ascending_node,
            self.mean_anomaly,
            self.orbital_period,
        )
    }
}

pub fn calculate_keplerian_orbital_elements(
    gm: f64,
    position: Axes,
    velocity: Axes,
) -> KeplerianElements {
    // Based on the implementation of Chambers in Mercury
    // Calculates Keplerian orbital elements given relative coordinates and
    // velocities, and GM = G times the sum of the masses.

    let (x, y, z) = position.unpack();
    let (u, v, w) = velocity.unpack();

    let mut eccentricity;
    let mut inclination = 0.;

    //  Local
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

    let semi_major_axis = gm * r / (2.0 * gm - r * v2);

    // Inclination and longitude of ascending node
    let mut longitude_of_ascending_node;
    let cos_i = (hz / h).clamp(-1., 1.);
    if abs!(cos_i) < 1. {
        inclination = acos!(cos_i);
        if hy > 0. {
            inclination = TWO_PI - inclination;
        }
        longitude_of_ascending_node = hx.atan2(-hy);
        if longitude_of_ascending_node < 0. {
            longitude_of_ascending_node += TWO_PI;
        }
        if inclination == 0. {
            longitude_of_ascending_node = 0.;
        }
    } else {
        if cos_i > 0. {
            inclination = 0.;
        }
        if cos_i < 0. {
            inclination = PI;
        }
        longitude_of_ascending_node = 0.;
    }

    // Eccentricity and perihelion distance
    let temp = 1. + s * (v2 / gm - 2. / r);
    if temp <= 0. {
        eccentricity = 0.;
    } else {
        eccentricity = sqrt!(temp);
    }
    let perihelion_distance = s / (1. + eccentricity);
    if eccentricity < 3.0e-8 {
        eccentricity = 0.;
    }

    let argument_perihelion = if eccentricity == 0. {
        0.
    } else {
        let ex_rot =
            ex * cos!(longitude_of_ascending_node) + ey * sin!(longitude_of_ascending_node);
        let ey_rot =
            -ex * sin!(longitude_of_ascending_node) + ey * cos!(longitude_of_ascending_node);
        let ez_rot = ez;
        let ex_rot_2 = ex_rot;
        let ey_rot_2 = ey_rot * cos!(inclination) + ez_rot * sin!(inclination);

        atan2!(ey_rot_2, ex_rot_2)
    };

    let longitude_perihelion = argument_perihelion + longitude_of_ascending_node;

    // True anomaly
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

    // Mean anomaly
    let mut mean_anomaly;
    if eccentricity == 0. {
        mean_anomaly = true_anomaly;
    } else {
        let mut cos_bige = ((1. / eccentricity) * (1. - (r / semi_major_axis))).clamp(-1., 1.);

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
        mean_anomaly = mean_anomaly.rem_euclid(TWO_PI);
    }

    // Given relative coordinates and velocities (of the body 1 respect to body 2),
    // and GM = G times the sum of the masses (body 1 + body 2)
    let orbital_period = (TWO_PI / sqrt!(gm)) * semi_major_axis.powf(3. / 2.);

    KeplerianElements {
        semi_major_axis,
        perihelion_distance,
        eccentricity,
        inclination,
        longitude_perihelion,
        longitude_of_ascending_node,
        mean_anomaly,
        orbital_period,
    }
}

// Based on the implementation of Chambers in Mercury
// Calculates Keplerian orbital elements given relative coordinates and
// velocities, and GM = G times the sum of the masses.
pub fn calculate_perihelion_distance_and_eccentricity(
    gm: f64,
    position: Axes,
    velocity: Axes,
) -> (f64, f64) {
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
    let perihelion_distance = s / (1. + eccentricity);

    (perihelion_distance, eccentricity)
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
        // Transform target to percent as in transforming x[left]..x[right] to 0..1
        let target_percent = (target_x - x_left) / (x[right] - x_left);
        y[left] * (1. - target_percent) + y[right] * target_percent
    };

    (target_y, left)
}

pub fn find_indices_around_target_value(data: &[f64], target_value: f64) -> (usize, usize) {
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
