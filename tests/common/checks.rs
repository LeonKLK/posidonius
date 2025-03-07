use posidonius::{
    EvolutionType, Particle, RotationalFlatteningEffect, RotationalFlatteningModel, TidalModel,
    TidesEffect,
};
use time::{OffsetDateTime, format_description};

// If creep coplanar tides and rotational flattening are set, both
// need to have the same uniform viscosity coefficient parameter
pub fn check_uniform_viscosity_coefficient(particle: &Particle) {
    let disabled_tides = matches!(particle.tides.effect, TidesEffect::Disabled);
    let disabled_rotational_flattening = matches!(
        particle.rotational_flattening.effect,
        RotationalFlatteningEffect::Disabled
    );
    if !disabled_tides && !disabled_rotational_flattening {
        let (creep_coplanar_tides, particle_uniform_viscosity_coefficient_for_tides) = {
            if let TidesEffect::CentralBody(TidalModel::CreepCoplanar(params))
            | TidesEffect::OrbitingBody(TidalModel::CreepCoplanar(params)) =
                &particle.tides.effect
            {
                (true, params.uniform_viscosity_coefficient)
            } else {
                (false, 0.)
            }
        };
        let (
            creep_coplanar_rotational_flattening,
            particle_uniform_viscosity_coefficient_for_rotational_flattenning,
        ) = {
            if let RotationalFlatteningEffect::CentralBody(
                RotationalFlatteningModel::CreepCoplanar(params),
            )
            | RotationalFlatteningEffect::OrbitingBody(
                RotationalFlatteningModel::CreepCoplanar(params),
            ) = &particle.rotational_flattening.effect
            {
                (true, params.uniform_viscosity_coefficient)
            } else {
                (false, 0.)
            }
        };
        if (creep_coplanar_tides && !creep_coplanar_rotational_flattening)
            || (!creep_coplanar_tides && creep_coplanar_rotational_flattening)
        {
            panic!(
                "[ERROR {} UTC] When using Creep Coplanar Tidal or rotational flattening effects, both effects need to be Creep Coplanar and not just one of them (e.g., it cannot be mixed with ConstantTimeLag or OblateSpheroid).",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        } else if creep_coplanar_tides && creep_coplanar_rotational_flattening {
            let diff_uniform_viscosity_coefficient = f64::abs(
                particle_uniform_viscosity_coefficient_for_tides
                    - particle_uniform_viscosity_coefficient_for_rotational_flattenning,
            );
            assert!(
                diff_uniform_viscosity_coefficient <= 1.0e-16,
                "[ERROR {} UTC] When using Creep Coplanar Tidal and rotational flattening effects, the uniform viscosity coefficient must be identical {:.16}.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                diff_uniform_viscosity_coefficient
            );
        }
    }
}

pub fn evolution_warnings(evolution: EvolutionType) {
    match evolution {
        EvolutionType::GalletBolmont2017(_) => {
            println!(
                "[WARNING {} UTC] Bodies with GalletBolmont2017 evolution will ignore initial radius and dissipation factor.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            println!(
                "[WARNING {} UTC] GalletBolmont2017 prescription theoretically only works for circular orbits and non inclined orbits, use carefully.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        EvolutionType::BolmontMathis2016(_) => {
            println!(
                "[WARNING {} UTC] Bodies with Baraffe2015 evolution will ignore initial radius and radius of gyration.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            println!(
                "[WARNING {} UTC] BolmontMathis2016 prescription theoretically only works for circular orbits and non inclined orbits, use carefully. ",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        EvolutionType::Baraffe2015(_) => println!(
            "[WARNING {} UTC]  ",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        ),
        EvolutionType::Leconte2011(_) => println!(
            "[WARNING {} UTC] Bodies with Leconte2011 evolution will ignore initial radius and radius of gyration.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        ),
        EvolutionType::Baraffe1998(_) => println!(
            "[WARNING {} UTC] Bodies with Baraffe1998 evolution will ignore initial radius. ",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        ),
        EvolutionType::LeconteChabrier2013(false) => println!(
            "[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration and love number.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        ),
        EvolutionType::LeconteChabrier2013(true) => {
            println!(
                "[WARNING {} UTC] Bodies with Jupiter evolution will ignore initial radius, radius of gyration, love number and dissipation factor.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            println!(
                "[WARNING {} UTC] LeconteChabrier2013(true) prescription theoretically only works for circular orbits and non inclined orbits, use carefully.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
        }
        EvolutionType::NonEvolving => {}
    }
}
