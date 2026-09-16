// Implemented by Alexandre Revol alexandre.revol@unige.ch
use super::TidesEffect;
use crate::constants::{DAY, G, TWO_PI};
use crate::tools;
use crate::tools::KeplerianElements;
use crate::{Axes, Particle};
use serde::{Deserialize, Serialize};

mod polynomials;
pub use polynomials::Polynomials;

mod love_number;
pub use love_number::LoveNumber;

mod components_2d;
use components_2d::calculate_2d_tidal_force_components;

mod components_3d;
use components_3d::calculate_3d_tidal_force_components;

// Structure declare the intput parameters needed for the Kaula model calculation.
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct KaulaParameters {
    // Externally provided spectrum for excitation frequency, real and imaginary part.
    #[serde(default)]
    pub love_numbers: LoveNumber,
    #[serde(default)]
    pub tidal_force: Axes,

    // Internally calculated kaula polynomials for mpq.
    #[serde(skip)]
    pub polynomials: Polynomials,
}

pub fn calculate_tidal_force(
    tidal_host_particle: &mut Particle,
    particle: &mut Particle,
    central_body: bool,
) -> Axes {
    //// TODO: Reconsider how to compute the stellar tide here instead of in tides.rs while allowing a mixture of tidal models
    let gm = G * (tidal_host_particle.mass + particle.mass);

    let (obliquity, keplerian_elements) = if central_body {
        (
            tools::calculate_inclination_orbital_equatorial_plane(
                tidal_host_particle.heliocentric_position,
                tidal_host_particle.heliocentric_velocity,
                particle.spin,
            ),
            tools::calculate_keplerian_orbital_elements(
                gm,
                tidal_host_particle.heliocentric_position,
                tidal_host_particle.heliocentric_velocity,
            ),
        )
    } else {
        (
            tools::calculate_inclination_orbital_equatorial_plane(
                particle.heliocentric_position,
                particle.heliocentric_velocity,
                particle.spin,
            ),
            tools::calculate_keplerian_orbital_elements(
                gm,
                particle.heliocentric_position,
                particle.heliocentric_velocity,
            ),
        )
    };

    let (tidal_force, tidal_force_secular) = calculate_2d_or_3d_tidal_force_components(
        tidal_host_particle,
        particle,
        central_body,
        obliquity,
        &keplerian_elements,
    );

    calculate_tidal_force_component(
        tidal_host_particle,
        particle,
        central_body,
        tidal_force,
        tidal_force_secular,
    )
}

fn calculate_tidal_force_component(
    tidal_host_particle: &mut Particle,
    particle: &mut Particle,
    central_body: bool,
    tidal_force: (f64, f64, f64),
    tidal_force_secular: (f64, f64, f64),
) -> Axes {
    // Spherical component of the tidal force
    // The radial component is the force applicated through the radial axis
    // The normal component act on the co longitude axis
    // The orthogonal component act on the co latitude axis
    // Required for denergy_dt calculation

    let (normal_component, orthogonal_component, radial_component) = tidal_force;

    let (normal_component_secular, orthogonal_component_secular, radial_component_secular) =
        tidal_force_secular;

    if central_body {
        tidal_host_particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide = orthogonal_component;
    } else {
        particle
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide = orthogonal_component;
    }

    // Cartesian tidal force computed by projection of the spherical coordinates
    let components = cartesian_projection_of_spherical_coordinates(
        tidal_host_particle,
        particle,
        central_body,
        normal_component,
        orthogonal_component,
        radial_component,
    );

    // Secular part of the tidal torque (simplified from the rapid varying phases)
    let secular_projection = cartesian_projection_of_spherical_coordinates(
        tidal_host_particle,
        particle,
        central_body,
        normal_component_secular,
        orthogonal_component_secular,
        radial_component_secular,
    );

    if central_body {
        // Stellar tide: `particle` is the star, whose single KaulaParameters::tidal_force is
        // shared by all planets and would only keep the last planet's force. Store the secular
        // force of this star-planet pair on the planet (`tidal_host_particle`) instead.
        tidal_host_particle
            .tides
            .parameters
            .internal
            .kaula_stellar_tide_secular_force = secular_projection;
    } else {
        particle.tides.get_kaula_mut().tidal_force = secular_projection;
    }

    components
}

fn calculate_2d_or_3d_tidal_force_components(
    tidal_host_particle: &Particle,
    particle: &mut Particle,
    central_body: bool,
    obliquity: f64,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    // !central_body ==> planetary tide ==> particle is the planet
    // central_body ==> stellar tide ==> tidal_host_particle is the planet
    // Keplerian elements
    let eccentricity = keplerian_elements.eccentricity;
    let orbital_period = keplerian_elements.orbital_period;

    // Planetary spin in [rad.s^-1]
    let spin = sqrt!(particle.norm_spin_vector_2) / DAY;
    // Orbital mean motion in [rad.s^-1]
    let orbital_frequency = TWO_PI / (orbital_period * DAY);

    let kaula = particle.tides.get_kaula_mut();
    // Update the kaula eccentricity, used by both 2D and 3D
    kaula.polynomials.update_eccentricity_2d(eccentricity);

    // Update the love number values and compute the tidal force components.
    // Decide if this is a 2d or 3d case
    if obliquity <= 1.0e-8 {
        // 2D case
        // refresh only love numbers used by 2D calculations
        kaula.love_numbers.refresh_cache_partial(
            central_body,
            spin,
            orbital_frequency,
            eccentricity,
        );

        calculate_2d_tidal_force_components(
            tidal_host_particle,
            particle,
            central_body,
            keplerian_elements,
        )
    } else {
        // 3D case
        kaula
            .love_numbers
            .refresh_cache_full(central_body, spin, orbital_frequency, eccentricity);

        // Update the kaula inclination and 3D eccentricty, used only by 3D
        kaula.polynomials.update_inclination(obliquity);
        kaula.polynomials.update_eccentricity_3d(eccentricity);

        calculate_3d_tidal_force_components(
            tidal_host_particle,
            particle,
            central_body,
            keplerian_elements,
        )
    }
}

// TODO physicist rename
fn trig_angles(particle: &Particle) -> (f64, f64, f64, f64) {
    let radial_distance = particle.tides.parameters.internal.distance;
    let (x, y, z) = particle.tides.coordinates.position.unpack();

    let coplanar_distance = sqrt!(x.powi(2) + y.powi(2));
    let cos_phi = x / coplanar_distance;
    let sin_phi = y / coplanar_distance;
    let cos_theta = z / radial_distance;
    let sin_theta = coplanar_distance / radial_distance;

    (cos_phi, cos_theta, sin_phi, sin_theta)
}

fn cartesian_projection_of_spherical_coordinates(
    tidal_host_particle: &Particle,
    particle: &mut Particle,
    central_body: bool,
    normal_component: f64,
    orthogonal_component: f64,
    radial_component: f64,
) -> Axes {
    // Spherical coordinate
    // The following elements correspond to the coordinate in the spherical coordinate
    // The coplanar distance is the radial distance projected in the x-y plane
    // The theta angle is the angle between the radial distance vector with respect to the z axis
    // The phi angle is the angle between the coplanar distance with respect to the x axis
    let (cos_phi, cos_theta, sin_phi, sin_theta) = if central_body {
        trig_angles(tidal_host_particle)
    } else {
        trig_angles(particle)
    };

    Axes::from(
        radial_component * sin_theta * cos_phi
        + normal_component * cos_theta * cos_phi
        - orthogonal_component * sin_phi,

        radial_component * sin_theta * sin_phi
        + normal_component * cos_theta * sin_phi
        + orthogonal_component * cos_phi,

        radial_component * cos_theta
        - normal_component * sin_theta
    )
}

// Calculate tidal torque due to tidal forces
pub fn calculate_torque_due_to_tides(
    tidal_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    let mut position = particle.tides.coordinates.position;

    let tidal_force = if central_body {
        if matches!(
            &tidal_host_particle.tides.effect,
            TidesEffect::CentralBody(_)
        ) {
            // additive inversion of positions if the host particle is the central body
            position.negate();
        }
        // If this is the central body, take the secular kaula force of this star-planet pair,
        // which calculate_tidal_force_component stored on the planet (`particle`)
        particle
            .tides
            .parameters
            .internal
            .kaula_stellar_tide_secular_force
    } else {
        // If it is not the central body, take the kaula tidal froce from the other particle
        particle.tides.get_kaula().tidal_force
    };
    // The negative sign of the position vector for stellar tide: Torque = r cross F
    // The r here should be the vector point from th primary (perturber) to the secondary (perturbed)
    // Thus we added a minus sign here as particle.tides.coordinates.position Axes are always heliocentric.

    // Let the torque be the cross product of the radial distance vector and the tidal force vector
    let torque_due_to_tides_x = position.y * tidal_force.z - position.z * tidal_force.y;
    let torque_due_to_tides_y = position.z * tidal_force.x - position.x * tidal_force.z;
    let torque_due_to_tides_z = position.x * tidal_force.y - position.y * tidal_force.x;

    Axes::from(
        torque_due_to_tides_x,
        torque_due_to_tides_y,
        torque_due_to_tides_z,
    )
}

pub fn select_eccentricty_order_q(eccentricity: f64) -> (usize, usize) {
    match () {
        // Select the order of the summation q over the eccentricity function G_lpq
        () if (eccentricity > 0.30) => (0, 15), // q: -7 <= q <= 7
        () if (eccentricity > 0.25) => (1, 14), // q: -6 <= q <= 6
        () if (eccentricity > 0.20) => (2, 13), // q: -5 <= q <= 5
        () if (eccentricity > 0.15) => (3, 12), // q: -4 <= q <= 4
        () if (eccentricity > 0.10) => (4, 11), // q: -3 <= q <= 3
        () if (eccentricity <= 0.1) => (5, 10), // q: -2 <= q <= 2
        () => unreachable!(),
    }
}

// Calculate the phases of the 2-Kaula transformed tidal forces
fn alpha_pqkj(
    p: usize,
    q: usize,
    k: usize,
    j: usize,
    mean_anomaly: f64,
    argument_perihelion: f64,
) -> f64 {
    let (p, q, k, j) = (p as i32, q as i32, k as i32, j as i32);
    f64!(2 * p - 2 * k + j - q) * mean_anomaly + f64!(2 * (p - k)) * argument_perihelion
}

/// Prefactor of the kaula tidal force: G m_perturber^2 R_deformed^5 / (a^6 r).
/// `tidal_deformed_body` is the body raising the tidal bulge (its radius enters), and
/// `tidal_perturber` is the body raising the tide (its mass enters):
/// - planetary tide: deformed body = planet, perturber = star;
/// - stellar tide:   deformed body = star,   perturber = planet.
fn calculate_base_constant(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    heliocentric_radius: f64,
    semi_major_axis: f64,
) -> f64 {
    (G * tidal_perturber.mass.powi(2) * tidal_deformed_body.radius.powi(5))
        / (semi_major_axis.powi(6) * heliocentric_radius)
}

/// 2D prefactor: `calculate_base_constant` divided by sin(theta) of the planet on its orbit.
/// `orbit_position` must be the planet's position (tides coordinates) for both the planetary
/// and the stellar tide, since the angles always describe the planet's orbit.
fn calculate_2d_constant(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    orbit_position: Axes,
    heliocentric_radius: f64,
    semi_major_axis: f64,
) -> f64 {
    let distance = orbit_position.norm();
    let (x, y, _z) = orbit_position.unpack();
    let coplanar_distance = sqrt!(x.powi(2) + y.powi(2));
    let sin_theta = coplanar_distance / distance;

    -(calculate_base_constant(
        tidal_deformed_body,
        tidal_perturber,
        heliocentric_radius,
        semi_major_axis,
    ) / sin_theta)
}
