use super::Axes;
use crate::constants::K2;
use crate::effects::{
    Disk, DiskEffect, EvolutionType, GeneralRelativity, GeneralRelativityEffect,
    RotationalFlattening, RotationalFlatteningEffect, Tides, TidesEffect, Wind, WindEffect,
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum Reference {
    MostMassiveParticle,
    Particle(usize), // Index
}

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub struct Particle {
    pub id: usize, // Unique internal identifier, to be set by the universe
    pub mass: f64,
    pub mass_g: f64,
    pub radius: f64,
    // Inertial frame where the center of mass of the system is at rest with respect to the origin of the coordinate system
    // (i.e., barycentric frame)
    pub inertial_position: Axes,
    pub inertial_velocity: Axes,
    pub inertial_acceleration: Axes,
    pub inertial_additional_acceleration: Axes,
    // Positions/velocities in a heliocentric frame: the host is at rest with respect to the origin of the coordinate system
    // where the host is the most massive particle in the universe
    pub heliocentric_position: Axes,
    pub heliocentric_velocity: Axes,
    pub heliocentric_distance: f64,               // Optimization
    pub heliocentric_radial_velocity: f64,        // Optimization
    pub heliocentric_norm_velocity_vector: f64,   // Optimization
    pub heliocentric_norm_velocity_vector_2: f64, // Optimization
    // Spin
    pub spin: Axes,
    pub norm_spin_vector_2: f64,
    pub angular_momentum: Axes,
    pub dangular_momentum_dt: Axes, // Force
    pub radius_of_gyration_2: f64, // radius of gyration square can be computed in terms of the mass moment of inertia, which
    // depends on the shape of the body and determines the torque needed for a desired angular acceleration
    pub moment_of_inertia: f64, // Spin related
    pub reference: Reference,   // Particle of reference for computing keplerian orbital parameters
    pub tides: Tides,
    pub rotational_flattening: RotationalFlattening,
    pub general_relativity: GeneralRelativity,
    pub wind: Wind,
    pub disk: Disk,
    pub evolution: EvolutionType,
}

impl Default for Particle {
    fn default() -> Self {
        let k_factor = 0.;
        let rotation_saturation = 0.;
        Self {
            id: 0,
            mass: 0.,
            mass_g: 0.,
            radius: 0.,
            inertial_position: Axes::new(),
            inertial_velocity: Axes::new(),
            inertial_acceleration: Axes::new(),
            inertial_additional_acceleration: Axes::new(),
            heliocentric_position: Axes::new(),
            heliocentric_velocity: Axes::new(),
            heliocentric_distance: 0.,
            heliocentric_radial_velocity: 0.,
            heliocentric_norm_velocity_vector: 0.,
            heliocentric_norm_velocity_vector_2: 0.,
            spin: Axes::new(),
            norm_spin_vector_2: 0.,
            angular_momentum: Axes::new(),
            dangular_momentum_dt: Axes::new(),
            radius_of_gyration_2: 0.,
            moment_of_inertia: 0.,
            reference: Reference::MostMassiveParticle,
            tides: Tides::new(TidesEffect::Disabled),
            rotational_flattening: RotationalFlattening::new(RotationalFlatteningEffect::Disabled),
            general_relativity: GeneralRelativity::new(GeneralRelativityEffect::Disabled),
            wind: Wind::new(WindEffect::Disabled, k_factor, rotation_saturation),
            disk: Disk::new(DiskEffect::Disabled),
            evolution: EvolutionType::NonEvolving,
        }
    }
}

impl Particle {
    pub fn new(
        mass: f64,
        radius: f64,
        radius_of_gyration: f64,
        position: Axes,
        velocity: Axes,
        spin: Axes,
    ) -> Particle {
        // Default effects: None
        let radius_of_gyration_2 = radius_of_gyration.powi(2);
        let moment_of_inertia = mass * radius_of_gyration_2 * radius.powi(2);
        Particle {
            mass,
            mass_g: mass * K2,
            radius,
            heliocentric_position: position,
            heliocentric_velocity: velocity,
            spin,
            norm_spin_vector_2: (spin.x.powi(2)) + (spin.y.powi(2)) + (spin.z.powi(2)),
            angular_momentum: Axes::from(
                moment_of_inertia * spin.x,
                moment_of_inertia * spin.y,
                moment_of_inertia * spin.z,
            ),
            radius_of_gyration_2,
            moment_of_inertia,
            ..Default::default()
        }
    }

    pub fn update_spin(&mut self) {
        self.spin.x = self.angular_momentum.x / self.moment_of_inertia;
        self.spin.y = self.angular_momentum.y / self.moment_of_inertia;
        self.spin.z = self.angular_momentum.z / self.moment_of_inertia;
        // norm needed for rotational flattening (torque and accelerations) and evolution
        self.norm_spin_vector_2 =
            (self.spin.x.powi(2)) + (self.spin.y.powi(2)) + (self.spin.z.powi(2));
    }

    pub fn view(&self) -> ParticleView {
        ParticleView {
            angular_momentum: self.angular_momentum,
            inertial_velocity: self.inertial_velocity,
        }
    }
}

// Instead of cloning the full particle, only clone the used fields
#[derive(Debug, Copy, Clone)]
pub struct ParticleView {
    pub angular_momentum: Axes,
    pub inertial_velocity: Axes,
}
