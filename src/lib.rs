#[macro_use]
extern crate math_macros;

pub mod constants;
mod effects;
mod integrator;
mod particles;
pub mod tools;

pub use effects::tides::kaula::LoveNumber;
pub use effects::{
    ConstantTimeLagParameters, CreepCoplanarParameters, Disk, DiskEffect, DiskProperties,
    EvolutionType, Evolver, GeneralRelativity, GeneralRelativityEffect,
    GeneralRelativityImplementation, KaulaParameters, OblateSpheroidParameters, Polynomials,
    RotationalFlattening, RotationalFlatteningEffect, RotationalFlatteningModel, TidalModel,
    TideComposition, Tides, TidesEffect, Wind, WindEffect,
};
pub use integrator::*;
pub use particles::{Axes, ConsiderEffects, IgnoreGravityTerms, Particle, Reference, Universe};
