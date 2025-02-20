use posidonius;

// NewFeatures:solar_like_for_kaula_1, solar_like_for_kaula_2, red_dwarf_like_for_kaula
#[allow(dead_code)]
pub fn solar_like_for_kaula_1(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    _general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::BolmontMathis2016(_) => {}
        posidonius::EvolutionType::GalletBolmont2017(_) => {}
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"
            );
        }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration = 0.2645751311; // Sun
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period = 1.; // days
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::DisabledModel,
    ));

    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    let star_general_relativity =
        posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::Disabled);
    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Disabled,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn solar_like_for_kaula_2(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    _general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::BolmontMathis2016(_) => {}
        posidonius::EvolutionType::GalletBolmont2017(_) => {}
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"
            );
        }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration = 0.2645751311; // Sun
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period = 2.9; // days
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::DisabledModel,
    ));

    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    let star_general_relativity =
        posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::Disabled);
    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Disabled,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn red_dwarf_like_for_kaula(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    _general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::BolmontMathis2016(_) => {}
        posidonius::EvolutionType::GalletBolmont2017(_) => {}
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"
            );
        }
    }
    let radius_factor: f64 = 0.117;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration = 4.47e-01; // Brown dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period = 1.; // days
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::DisabledModel,
    ));

    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    let star_general_relativity =
        posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::Disabled);
    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Disabled,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn solar_like(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::BolmontMathis2016(_) => {}
        posidonius::EvolutionType::GalletBolmont2017(_) => {}
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"
            );
        }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration: f64 = 2.43e-01; // Sun
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period: f64 = 8.; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period / 24.); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_love_number: f64 = 0.03; // Sun
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 4.992 * 3.845764e-2; // -66+64
    let star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: star_dissipation_factor,
        dissipation_factor_scale: star_dissipation_factor_scale,
        love_number: star_love_number,
    };
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::ConstantTimeLag(star_tidal_model_params),
    ));
    let star_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters {
        love_number: star_love_number,
    };
    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::CentralBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                star_rotational_flattening_oblate_spheroid_params,
            ),
        ));
    let star_general_relativity = posidonius::GeneralRelativity::new(
        posidonius::GeneralRelativityEffect::CentralBody(general_relativity_implementation),
    );
    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let star_wind_k_factor = 4.0e-18;
    let star_wind_rotation_saturation = 1.7592918860102842;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Interaction,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn solar_like_with_disk(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::BolmontMathis2016(_) => {}
        posidonius::EvolutionType::GalletBolmont2017(_) => {}
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015, BolmontMathis2016, GalletBolmont2017 or non evolving to create a solar like body!"
            );
        }
    }
    let radius_factor: f64 = 1.0;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration: f64 = 2.43e-01; // Sun
    //
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period: f64 = 8.; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period / 24.); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_love_number: f64 = 0.03; // Sun
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 4.992 * 3.845764e-2; // -66+64
    let star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: star_dissipation_factor,
        dissipation_factor_scale: star_dissipation_factor_scale,
        love_number: star_love_number,
    };
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::ConstantTimeLag(star_tidal_model_params),
    ));
    let star_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters {
        love_number: star_love_number,
    };
    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::CentralBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                star_rotational_flattening_oblate_spheroid_params,
            ),
        ));
    let star_general_relativity = posidonius::GeneralRelativity::new(
        posidonius::GeneralRelativityEffect::CentralBody(general_relativity_implementation),
    );
    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let star_wind_k_factor = 4.0e-18;
    let star_wind_rotation_saturation = 1.7592918860102842;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Interaction,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    // Disk
    let surface_density_normalization_gcm = 1000.; // g.cm^-2
    let surface_density_normalization_si = surface_density_normalization_gcm * 1.0e-3 * 1.0e4; // kg.m^-2
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::CentralBody(
        posidonius::DiskProperties {
            inner_edge_distance: 0.01,  // AU
            outer_edge_distance: 100.0, // AU
            lifetime: 1.0e5 * 365.25e0, // days
            alpha: 1.0e-2,
            surface_density_normalization: surface_density_normalization_si
                * (1.0 / posidonius::constants::M_SUN)
                * posidonius::constants::AU.powi(2), // Msun.AU^-2
            mean_molecular_weight: 2.4,
        },
    ));
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn brown_dwarf(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::Leconte2011(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Leconte2011 or non evolving to create a Brown Dwarf body!"
            );
        }
    }

    let (star_rotation_period, star_love_number) = match star_evolution {
        posidonius::EvolutionType::NonEvolving => {
            let rotation_period: f64 = 70.0; // hours
            let love_number: f64 = 0.307; // BrownDwarf
            (rotation_period, love_number)
        }
        posidonius::EvolutionType::Leconte2011(mass) => {
            let rotation_period; // hours
            let love_number;
            if (0.0099..=0.0101).contains(&mass) {
                rotation_period = 8.0;
                love_number = 0.3790;
            } else if (0.0119..=0.0121).contains(&mass) {
                rotation_period = 13.0;
                love_number = 0.3780;
            } else if (0.0149..=0.0151).contains(&mass) {
                rotation_period = 19.0;
                love_number = 0.3760;
            } else if (0.0199..=0.0201).contains(&mass) {
                rotation_period = 24.0;
                love_number = 0.3690;
            } else if (0.0299..=0.0301).contains(&mass) {
                rotation_period = 30.0;
                love_number = 0.3550;
            } else if (0.0399..=0.0401).contains(&mass) {
                rotation_period = 36.0;
                love_number = 0.3420;
            } else if (0.0499..=0.0501).contains(&mass) {
                rotation_period = 41.0;
                love_number = 0.3330;
            } else if (0.0599..=0.0601).contains(&mass) {
                rotation_period = 47.0;
                love_number = 0.3250;
            } else if (0.0699..=0.0701).contains(&mass) {
                rotation_period = 53.0;
                love_number = 0.3110;
            } else if (0.0719..=0.0721).contains(&mass) {
                rotation_period = 58.0;
                love_number = 0.3080;
            } else if (0.0749..=0.0751).contains(&mass) {
                rotation_period = 64.0;
                love_number = 0.3070;
            } else if (0.0799..=0.0801).contains(&mass) {
                rotation_period = 70.0;
                love_number = 0.3070;
            } else {
                panic!("The evolution type Leconte2011 does not support a mass of {mass} Msun!");
            }
            (rotation_period, love_number)
        }
        _ => {
            panic!(
                "Evolution type should be solar like or non evolving to create a solar like body!"
            );
        }
    };

    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006 * 3.845764e4; // -60+64
    // Radius of gyration
    let star_radius_of_gyration: f64 = 4.41e-01; // Brown dwarf
    //
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period / 24.); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: star_dissipation_factor,
        dissipation_factor_scale: star_dissipation_factor_scale,
        love_number: star_love_number,
    };
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::ConstantTimeLag(star_tidal_model_params),
    ));
    let star_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters {
        love_number: star_love_number,
    };
    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::CentralBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                star_rotational_flattening_oblate_spheroid_params,
            ),
        ));
    let star_general_relativity = posidonius::GeneralRelativity::new(
        posidonius::GeneralRelativityEffect::CentralBody(general_relativity_implementation),
    );
    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Disabled,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn m_dwarf(
    star_mass: f64,
    star_evolution: posidonius::EvolutionType,
    general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    match star_evolution {
        posidonius::EvolutionType::Baraffe1998(_) => {}
        posidonius::EvolutionType::Baraffe2015(_) => {}
        posidonius::EvolutionType::NonEvolving => {}
        _ => {
            panic!(
                "Evolution type should be Baraffe1998, Baraffe2015 or non evolving to create a M-Dwarf body!"
            );
        }
    }

    let radius_factor: f64 = 0.845649342247916;
    let star_radius: f64 = radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let star_radius_of_gyration: f64 = 4.47e-01; // M-dwarf
    // To calculate tidal forces, it is needed to have the central body/star at [0,0,0] and without velocity or acceleration (heliocentric)
    let star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let star_rotation_period: f64 = 70.0; // hours
    let star_angular_frequency = posidonius::constants::TWO_PI / (star_rotation_period / 24.); // days^-1
    let star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: star_angular_frequency,
    };

    let star_love_number: f64 = 0.307; // Brown Dwarf / M Dwarf
    // Disipation factor (sigma)
    let star_dissipation_factor_scale: f64 = 1.;
    // BD, Mdwarf: sigmast = 2.006d-60 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let star_dissipation_factor: f64 = 2.006 * 3.845764e4; // -60+64
    let star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: star_dissipation_factor,
        dissipation_factor_scale: star_dissipation_factor_scale,
        love_number: star_love_number,
    };
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::ConstantTimeLag(star_tidal_model_params),
    ));
    let star_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters {
        love_number: star_love_number,
    };
    let star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::CentralBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                star_rotational_flattening_oblate_spheroid_params,
            ),
        ));
    let star_general_relativity = posidonius::GeneralRelativity::new(
        posidonius::GeneralRelativityEffect::CentralBody(general_relativity_implementation),
    );
    let star_wind_k_factor = 0.;
    let star_wind_rotation_saturation = 0.;
    let star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Disabled,
        star_wind_k_factor,
        star_wind_rotation_saturation,
    );
    let star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut star = posidonius::Particle::new(
        star_mass,
        star_radius,
        star_radius_of_gyration,
        star_position,
        star_velocity,
        star_spin,
    );
    star.set_tides(star_tides);
    star.set_rotational_flattening(star_rotational_flattening);
    star.set_general_relativity(star_general_relativity);
    star.set_wind(star_wind);
    star.set_disk(star_disk);
    star.set_evolution(star_evolution);
    star
}

#[allow(dead_code)]
pub fn solar_like_primary(
    primary_star_mass: f64,
    primary_star_evolution: posidonius::EvolutionType,
) -> posidonius::Particle {
    let primary_star_radius_factor: f64 = 1.0;
    let primary_star_radius: f64 = primary_star_radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let primary_star_radius_of_gyration: f64 = 2.43e-01; // Sun
    //
    let primary_star_position = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    let primary_star_velocity = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: 0.,
    };
    // Initialization of stellar spin
    let primary_star_rotation_period: f64 = 8.; // hours
    let primary_star_angular_frequency =
        posidonius::constants::TWO_PI / (primary_star_rotation_period / 24.); // days^-1
    let primary_star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: primary_star_angular_frequency,
    };

    let primary_star_love_number: f64 = 0.03; // Sun
    // Disipation factor (sigma)
    let primary_star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-primary_star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let primary_star_dissipation_factor: f64 = 4.992 * 3.845764e-2; // -66+64
    let primary_star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: primary_star_dissipation_factor,
        dissipation_factor_scale: primary_star_dissipation_factor_scale,
        love_number: primary_star_love_number,
    };
    let primary_star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::ConstantTimeLag(primary_star_tidal_model_params),
    ));
    let primary_star_rotational_flattening_oblate_spheroid_params =
        posidonius::OblateSpheroidParameters {
            love_number: primary_star_love_number,
        };
    let primary_star_rotational_flattening =
        posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::CentralBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                primary_star_rotational_flattening_oblate_spheroid_params,
            ),
        ));
    let primary_star_general_relativity =
        posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::OrbitingBody);
    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let primary_star_wind_k_factor = 4.0e-18;
    let primary_star_wind_rotation_saturation = 1.7592918860102842;
    let primary_star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Interaction,
        primary_star_wind_k_factor,
        primary_star_wind_rotation_saturation,
    );
    let primary_star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut primary_star = posidonius::Particle::new(
        primary_star_mass,
        primary_star_radius,
        primary_star_radius_of_gyration,
        primary_star_position,
        primary_star_velocity,
        primary_star_spin,
    );
    primary_star.set_tides(primary_star_tides);
    primary_star.set_rotational_flattening(primary_star_rotational_flattening);
    primary_star.set_general_relativity(primary_star_general_relativity);
    primary_star.set_wind(primary_star_wind);
    primary_star.set_disk(primary_star_disk);
    primary_star.set_evolution(primary_star_evolution);
    primary_star
}

#[allow(dead_code)]
pub fn solar_like_secondary(
    primary_star: &posidonius::Particle,
    secondary_star_mass: f64,
    secondary_star_evolution: posidonius::EvolutionType,
    general_relativity_implementation: posidonius::GeneralRelativityImplementation,
) -> posidonius::Particle {
    let secondary_star_radius_factor: f64 = 1.0;
    let secondary_star_radius: f64 = secondary_star_radius_factor * posidonius::constants::R_SUN;
    // Radius of gyration
    let secondary_star_radius_of_gyration: f64 = 2.43e-01; // Sun

    let semimajor_axis = 100.;
    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    //let a: f64 = 0.018;                             // semi-major axis (in AU)
    let a: f64 = semimajor_axis; // semi-major axis (in AU)
    let e: f64 = 0.1; // eccentricity
    let i: f64 = 5. * posidonius::constants::DEG2RAD; // inclination (degrees)
    let mut p: f64 = 0. * posidonius::constants::DEG2RAD; // argument of pericentre (degrees)
    let n: f64 = 0. * posidonius::constants::DEG2RAD; // longitude of the ascending node (degrees)
    let l: f64 = 0. * posidonius::constants::DEG2RAD; // mean anomaly (degrees)
    p += n; // Convert to longitude of perihelion !!
    let q = a * (1.0 - e); // perihelion distance
    let gm: f64 = posidonius::constants::G * (primary_star.mass + secondary_star_mass);
    let (x, y, z, vx, vy, vz) =
        posidonius::tools::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let secondary_star_position = posidonius::Axes { x, y, z };
    let secondary_star_velocity = posidonius::Axes {
        x: vx,
        y: vy,
        z: vz,
    };

    // Initialization of stellar spin
    let secondary_star_rotation_period: f64 = 8.; // hours
    let secondary_star_angular_frequency =
        posidonius::constants::TWO_PI / (secondary_star_rotation_period / 24.); // days^-1
    let secondary_star_spin = posidonius::Axes {
        x: 0.,
        y: 0.,
        z: secondary_star_angular_frequency,
    };

    let secondary_star_love_number: f64 = 0.03; // Sun
    // Disipation factor (sigma)
    let secondary_star_dissipation_factor_scale: f64 = 1.;
    // Sun-like-secondary_star: sigmast = 4.992e-66 cgs, conversion to Msun-1.AU-2.day-1 = 3.845764022293d64
    let secondary_star_dissipation_factor: f64 = 4.992 * 3.845764e-2; // -66+64
    let secondary_star_tidal_model_params = posidonius::ConstantTimeLagParameters {
        dissipation_factor: secondary_star_dissipation_factor,
        dissipation_factor_scale: secondary_star_dissipation_factor_scale,
        love_number: secondary_star_love_number,
    };
    let secondary_star_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(
        posidonius::TidalModel::ConstantTimeLag(secondary_star_tidal_model_params),
    ));
    let secondary_star_rotational_flattening_oblate_spheroid_params =
        posidonius::OblateSpheroidParameters {
            love_number: secondary_star_love_number,
        };
    let secondary_star_rotational_flattening = posidonius::RotationalFlattening::new(
        posidonius::RotationalFlatteningEffect::OrbitingBody(
            posidonius::RotationalFlatteningModel::OblateSpheroid(
                secondary_star_rotational_flattening_oblate_spheroid_params,
            ),
        ),
    );
    let secondary_star_general_relativity = posidonius::GeneralRelativity::new(
        posidonius::GeneralRelativityEffect::CentralBody(general_relativity_implementation),
    );
    // Wind parametrisation (Bouvier 1997):
    //      wind_k_factor = 4.0e-18 # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    //      wind_rotation_saturation = 14. * TWO_PI/25.0 # = 1.7592918860102842, wsat in units of the spin of the Sun today
    let secondary_star_wind_k_factor = 4.0e-18;
    let secondary_star_wind_rotation_saturation = 1.7592918860102842;
    let secondary_star_wind = posidonius::Wind::new(
        posidonius::WindEffect::Interaction,
        secondary_star_wind_k_factor,
        secondary_star_wind_rotation_saturation,
    );
    let secondary_star_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut secondary_star = posidonius::Particle::new(
        secondary_star_mass,
        secondary_star_radius,
        secondary_star_radius_of_gyration,
        secondary_star_position,
        secondary_star_velocity,
        secondary_star_spin,
    );
    secondary_star.set_tides(secondary_star_tides);
    secondary_star.set_rotational_flattening(secondary_star_rotational_flattening);
    secondary_star.set_general_relativity(secondary_star_general_relativity);
    secondary_star.set_wind(secondary_star_wind);
    secondary_star.set_disk(secondary_star_disk);
    secondary_star.set_evolution(secondary_star_evolution);
    secondary_star
}
