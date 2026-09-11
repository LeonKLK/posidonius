// Implemented by Alexandre Revol alexandre.revol@unige.ch
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

// Roles in the kaula tidal force
// ------------------------------
// Every function below names its inputs by their physical role:
// - `tidal_deformed_body`: the body raising the tidal bulge that dissipates energy. Its spin,
//   radius and love number spectrum enter the force. Star for the stellar tide
//   (`central_body == true`), planet for the planetary tide (`central_body == false`).
// - `tidal_perturber`: the body whose gravity raises the bulge. Only its mass enters the force.
// - `orbit`: the star-planet orbit. It is ALWAYS taken from the planet (the orbiting body carries
//   the heliocentric position/velocity), whichever body is deformed.
// The `central_body` flag only selects the physics that genuinely differs between the two tides
// (love number parity, where the orthogonal component is stored, the sign conventions of the
// projection); it never decides which body provides the mass or the radius.

/// Orbital quantities of the star-planet pair, taken from the planet (the orbiting body).
#[derive(Debug, Copy, Clone)]
pub struct Orbit {
    /// Heliocentric position/velocity/distance of the planet (used for the keplerian elements).
    pub heliocentric_position: Axes,
    pub heliocentric_velocity: Axes,
    pub heliocentric_distance: f64,
    /// Position and distance of the planet as stored in `tides.coordinates` (used for the angles).
    pub position: Axes,
    pub distance: f64,
}

impl Orbit {
    pub fn from_planet(planet: &Particle) -> Self {
        Self {
            heliocentric_position: planet.heliocentric_position,
            heliocentric_velocity: planet.heliocentric_velocity,
            heliocentric_distance: planet.heliocentric_distance,
            position: planet.tides.coordinates.position,
            distance: planet.tides.parameters.internal.distance,
        }
    }
}

/// Tidal force on the tidally deformed body due to the perturber.
/// - Stellar tide (`central_body == true`): `tidal_deformed_body` is the star, `tidal_perturber`
///   the planet.
/// - Planetary tide (`central_body == false`): `tidal_deformed_body` is the planet,
///   `tidal_perturber` the star.
/// The orbit is always read from the planet. Both particles are mutable because the deformed
/// body caches its love numbers and secular force, and the planet stores the orthogonal component.
pub fn calculate_tidal_force(
    tidal_deformed_body: &mut Particle,
    tidal_perturber: &mut Particle,
    central_body: bool,
) -> Axes {
    //// TODO: Reconsider how to compute the stellar tide here instead of in tides.rs while allowing a mixture of tidal models
    let gm = G * (tidal_deformed_body.mass + tidal_perturber.mass);

    let orbit = Orbit::from_planet(planet(tidal_deformed_body, tidal_perturber, central_body));

    // Obliquity of the deformed body with respect to the orbit, and keplerian elements of the orbit.
    let obliquity = tools::calculate_inclination_orbital_equatorial_plane(
        orbit.heliocentric_position,
        orbit.heliocentric_velocity,
        tidal_deformed_body.spin,
    );
    let keplerian_elements = tools::calculate_keplerian_orbital_elements(
        gm,
        orbit.heliocentric_position,
        orbit.heliocentric_velocity,
    );

    let (tidal_force, tidal_force_secular) = calculate_2d_or_3d_tidal_force_components(
        tidal_deformed_body,
        tidal_perturber,
        &orbit,
        central_body,
        obliquity,
        &keplerian_elements,
    );

    calculate_tidal_force_component(
        tidal_deformed_body,
        tidal_perturber,
        &orbit,
        central_body,
        tidal_force,
        tidal_force_secular,
    )
}

/// The planet of the pair: the perturber for the stellar tide, the deformed body for the
/// planetary tide.
fn planet<'a>(
    tidal_deformed_body: &'a mut Particle,
    tidal_perturber: &'a mut Particle,
    central_body: bool,
) -> &'a mut Particle {
    if central_body {
        tidal_perturber
    } else {
        tidal_deformed_body
    }
}

fn calculate_tidal_force_component(
    tidal_deformed_body: &mut Particle,
    tidal_perturber: &mut Particle,
    orbit: &Orbit,
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

    // The orthogonal component is stored on the planet, in the slot of the tide that produced it.
    if central_body {
        tidal_perturber
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_stellar_tide = orthogonal_component;
    } else {
        tidal_deformed_body
            .tides
            .parameters
            .internal
            .orthogonal_component_of_the_tidal_force_due_to_planetary_tide = orthogonal_component;
    }

    // Cartesian tidal force computed by projection of the spherical coordinates
    let components = cartesian_projection_of_spherical_coordinates(
        orbit,
        normal_component,
        orthogonal_component,
        radial_component,
    );

    // Secular part of the tidal torque (simplified from the rapid varying phases)
    let secular_projection = cartesian_projection_of_spherical_coordinates(
        orbit,
        normal_component_secular,
        orthogonal_component_secular,
        radial_component_secular,
    );

    tidal_deformed_body.tides.get_kaula_mut().tidal_force = secular_projection;

    components
}

fn calculate_2d_or_3d_tidal_force_components(
    tidal_deformed_body: &mut Particle,
    tidal_perturber: &Particle,
    orbit: &Orbit,
    central_body: bool,
    obliquity: f64,
    keplerian_elements: &KeplerianElements,
) -> ((f64, f64, f64), (f64, f64, f64)) {
    // Keplerian elements
    let eccentricity = keplerian_elements.eccentricity;
    let orbital_period = keplerian_elements.orbital_period;

    // Spin of the tidally deformed body in [rad.s^-1]
    let spin = sqrt!(tidal_deformed_body.norm_spin_vector_2) / DAY;
    // Orbital mean motion in [rad.s^-1]
    let orbital_frequency = TWO_PI / (orbital_period * DAY);

    let kaula = tidal_deformed_body.tides.get_kaula_mut();
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
            tidal_deformed_body,
            tidal_perturber,
            orbit,
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
            tidal_deformed_body,
            tidal_perturber,
            orbit,
            central_body,
            keplerian_elements,
        )
    }
}

// TODO physicist rename
fn trig_angles(orbit: &Orbit) -> (f64, f64, f64, f64) {
    let radial_distance = orbit.distance;
    let (x, y, z) = orbit.position.unpack();

    let coplanar_distance = sqrt!(x.powi(2) + y.powi(2));
    let cos_phi = x / coplanar_distance;
    let sin_phi = y / coplanar_distance;
    let cos_theta = z / radial_distance;
    let sin_theta = coplanar_distance / radial_distance;

    (cos_phi, cos_theta, sin_phi, sin_theta)
}

fn cartesian_projection_of_spherical_coordinates(
    orbit: &Orbit,
    normal_component: f64,
    orthogonal_component: f64,
    radial_component: f64,
) -> Axes {
    // Spherical coordinate
    // The following elements correspond to the coordinate in the spherical coordinate
    // The coplanar distance is the radial distance projected in the x-y plane
    // The theta angle is the angle between the radial distance vector with respect to the z axis
    // The phi angle is the angle between the coplanar distance with respect to the x axis
    // The angles are those of the planet on its orbit, whichever body is deformed.
    let (cos_phi, cos_theta, sin_phi, sin_theta) = trig_angles(orbit);

    Axes::from(
        radial_component * sin_theta * cos_phi + normal_component * cos_theta * cos_phi
            - orthogonal_component * sin_phi,
        radial_component * sin_theta * sin_phi
            + normal_component * cos_theta * sin_phi
            + orthogonal_component * cos_phi,
        radial_component * cos_theta - normal_component * sin_theta,
    )
}

// Calculate tidal torque due to tidal forces
/// Torque on the tidally deformed body: r x F with F the secular kaula force cached on the
/// deformed body and r the position of the planet on its orbit.
pub fn calculate_torque_due_to_tides(
    tidal_deformed_body: &Particle,
    orbit: &Orbit,
    central_body: bool,
) -> Axes {
    let mut position = orbit.position;

    if central_body {
        // The negative sign of the position vector for stellar tide: Torque = r cross F
        // The r here should be the vector point from the primary (perturber) to the secondary (perturbed)
        // Thus we added a minus sign here as the orbit position is always heliocentric.
        position.negate();
    }
    let tidal_force = tidal_deformed_body.tides.get_kaula().tidal_force;

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

/// Prefactor of the kaula tidal force, `G m_perturber^2 R_deformed^5 / (a^6 r)`.
/// The radius is the one of the tidally deformed body, the mass the one of the perturber.
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
fn calculate_2d_constant(
    tidal_deformed_body: &Particle,
    tidal_perturber: &Particle,
    orbit: &Orbit,
    heliocentric_radius: f64,
    semi_major_axis: f64,
) -> f64 {
    let distance = orbit.position.norm();
    let (x, y, _z) = orbit.position.unpack();
    let coplanar_distance = sqrt!(x.powi(2) + y.powi(2));
    let sin_theta = coplanar_distance / distance;

    -(calculate_base_constant(
        tidal_deformed_body,
        tidal_perturber,
        heliocentric_radius,
        semi_major_axis,
    ) / sin_theta)
}
