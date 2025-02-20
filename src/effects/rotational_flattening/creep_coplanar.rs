// Implemented by Gabriel Oliveira Gomes
// Described in Gomes et al (2021): https://ui.adsabs.harvard.edu/abs/2021A%26A...651A..23G/abstract
use super::super::super::Axes;
use super::super::super::Particle;
use super::super::super::constants::K2;
use super::super::rotational_flattening::RotationalFlatteningEffect;
use super::super::rotational_flattening::RotationalFlatteningModel;
use super::super::tides::creep_coplanar::calculate_creep_coplanar_shape;

pub fn calculate_creep_coplanar_shapes(
    rotational_flattening_host_particle: &mut Particle,
    particles: &mut [Particle],
    more_particles: &mut [Particle],
) {
    for particle in particles.iter_mut().chain(more_particles.iter_mut()) {
        let particle_uniform_viscosity_coefficient = match particle.rotational_flattening.effect {
            RotationalFlatteningEffect::CentralBody(rotational_flattening_model)
            | RotationalFlatteningEffect::OrbitingBody(rotational_flattening_model) => {
                if let RotationalFlatteningModel::CreepCoplanar(params) =
                    rotational_flattening_model
                {
                    params.uniform_viscosity_coefficient
                } else {
                    0.
                }
            }
            RotationalFlatteningEffect::Disabled => 0.,
        };
        if particle_uniform_viscosity_coefficient == 0. {
            continue;
        }
        let consider_tides = false;
        let consider_rotational_flattening = true;
        let central_body = false;
        let shape = calculate_creep_coplanar_shape(
            rotational_flattening_host_particle,
            particle,
            consider_tides,
            consider_rotational_flattening,
            central_body,
        );
        particle.rotational_flattening.parameters.internal.shape = shape;
    }
}

pub fn calculate_torque_induced_by_rotational_flattening(
    rotational_flattening_host_particle: &Particle,
    particle: &Particle,
    central_body: bool,
) -> Axes {
    // Torque expression for studying creep tide tidal despinning (use torque.x = 0 instead to study stat. rotation, makes code run faster)
    let torque_induced_by_rotational_flattening_z = if central_body {
        // Host shape is planet dependent, it needs to be computed for each
        let consider_tides = false;
        let consider_rotational_flattening = true;
        let rotational_flattening_host_shape = calculate_creep_coplanar_shape(
            rotational_flattening_host_particle,
            particle,
            consider_tides,
            consider_rotational_flattening,
            central_body,
        );
        3.0 / 5.0
            * K2
            * rotational_flattening_host_particle.mass
            * particle.mass
            * rotational_flattening_host_particle.radius
            * rotational_flattening_host_particle.radius
            * rotational_flattening_host_shape.y
            / particle
                .rotational_flattening
                .parameters
                .internal
                .distance
                .powi(3)
    } else {
        // Planet shape is host dependent and it was already computed
        3.0 / 5.0
            * K2
            * particle.mass
            * rotational_flattening_host_particle.mass
            * particle.radius
            * particle.radius
            * particle.rotational_flattening.parameters.internal.shape.y
            / particle
                .rotational_flattening
                .parameters
                .internal
                .distance
                .powi(3)
    };

    Axes::from(0.0, 0.0, torque_induced_by_rotational_flattening_z)
}

pub fn calculate_acceleration_induced_by_rotational_flattering(
    rotational_flattening_host_particle: &Particle,
    particle: &Particle,
) -> Axes {
    let distance_4 = particle
        .rotational_flattening
        .parameters
        .internal
        .distance
        .powi(4);

    let factor = particle.rotational_flattening.parameters.internal.distance
        * (-3.0
            * K2
            * particle.mass
            * rotational_flattening_host_particle.mass
            * particle.radius.powi(2)
            / (5.0 * distance_4)
            * particle.rotational_flattening.parameters.internal.shape.z);

    // Expression for ONLY rotational flattenings
    let total_force_induced_by_rotational_flattening_x =
        particle.rotational_flattening.coordinates.position.x / factor;

    // Expression for ONLY rotational flattenings
    let total_force_induced_by_rotational_flattening_y =
        particle.rotational_flattening.coordinates.position.y / factor;

    Axes::from(
        total_force_induced_by_rotational_flattening_x,
        total_force_induced_by_rotational_flattening_y,
        0.0,
    )
}
