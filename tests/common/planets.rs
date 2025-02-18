extern crate posidonius;

// NewFeatures: new functions added: basic_planetary_configuration_for_stellar_kaula,
// basic_planetary_configuration_for_planetary_kaula, super_earth_like, trappist1_h_like
#[allow(dead_code)]
pub fn basic_planetary_configuration_for_stellar_kaula(star: &posidonius::Particle) -> Vec<posidonius::Particle> {
    let planet1_mass_factor: f64 = 5.0;
    let planet1_mass: f64 = planet1_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet1_evolution = posidonius::EvolutionType::NonEvolving;
    let planet1_semimajor_axis = 0.02;
    let planet1_eccentricity = 0.000001;
    let planet1_obliquity: f64 = 0. * posidonius::constants::DEG2RAD;
    let planet1 = super_earth_like(&star, planet1_mass, planet1_evolution, planet1_semimajor_axis, planet1_eccentricity, planet1_obliquity);

    return vec![planet1];
}

#[allow(dead_code)]
pub fn basic_planetary_configuration_for_planetary_kaula(star: &posidonius::Particle) -> Vec<posidonius::Particle> {
    let planet1_mass_factor: f64 = 0.3285168664;
    let planet1_mass: f64 = planet1_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet1_evolution = posidonius::EvolutionType::NonEvolving;
    let planet1_semimajor_axis = 0.0619224236258543;
    let planet1_eccentricity = 0.1;
    let planet1_obliquity: f64 = 0. * posidonius::constants::DEG2RAD;
    let planet1 = trappist1_h_like(&star, planet1_mass, planet1_evolution, planet1_semimajor_axis, planet1_eccentricity, planet1_obliquity);

    return vec![planet1];
}

#[allow(dead_code)]
pub fn basic_configuration(star: &posidonius::Particle) -> Vec<posidonius::Particle> {
    let planet1_mass_factor: f64 = 1.0;
    let planet1_mass: f64 = planet1_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet1_evolution = posidonius::EvolutionType::NonEvolving;
    let planet1_semimajor_axis = 0.5;
    let planet1_eccentricity = 0.1;
    let planet1_obliquity: f64 = 11.459156 * posidonius::constants::DEG2RAD; // 0.2 rad
    let planet1 = earth_like(&star, planet1_mass, planet1_evolution, planet1_semimajor_axis, planet1_eccentricity, planet1_obliquity);

    let planet2_mass_factor: f64 = 0.00095; // Earth masses
    let planet2_mass: f64 = planet2_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet2_evolution = posidonius::EvolutionType::LeconteChabrier2013(false); // Jupiter without dissipation of dynamical tides
    let planet2_semimajor_axis = 5.0;
    let planet2_eccentricity = 0.1;
    let planet2_obliquity: f64 = 45.0 * posidonius::constants::DEG2RAD;
    let planet2 = jupiter_like(&star, planet2_mass, planet2_evolution, planet2_semimajor_axis, planet2_eccentricity, planet2_obliquity);

    let planet3_mass_factor: f64 = 0.00095; // Earth masses
    let planet3_mass: f64 = planet3_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet3_evolution = posidonius::EvolutionType::NonEvolving;
    let planet3_semimajor_axis = 10.0;
    let planet3_eccentricity = 0.1; // this helps testing Kaula special case for small obliquities & eccentricity different than zero
    let planet3_obliquity: f64 = 0.0 * posidonius::constants::DEG2RAD; // 0.0 rad => this helps testing Kaula special case for small obliquities
    let planet3 = jupiter_like(&star, planet3_mass, planet3_evolution, planet3_semimajor_axis, planet3_eccentricity, planet3_obliquity);

    let planet4_mass_factor: f64 = 0.00095; // Earth masses
    let planet4_mass: f64 = planet4_mass_factor * posidonius::constants::M_EARTH; // Solar masses (3.0e-6 solar masses = 1 earth mass)
    let planet4_evolution = posidonius::EvolutionType::NonEvolving;
    let planet4_semimajor_axis = 15.0;
    let planet4_eccentricity = 0.0; // this helps testing Kaula special case for small obliquities & eccentricity exactly zero
    let planet4_obliquity: f64 = 0.0 * posidonius::constants::DEG2RAD; // 0.0 rad => this helps testing Kaula special case for small obliquities
    let planet4 = jupiter_like(&star, planet4_mass, planet4_evolution, planet4_semimajor_axis, planet4_eccentricity, planet4_obliquity);
    return vec![planet1, planet2, planet3, planet4];
}

#[allow(dead_code)]
pub fn super_earth_like(star: &posidonius::Particle, planet_mass: f64, planet_evolution: posidonius::EvolutionType, semimajor_axis: f64, eccentricity: f64, planet_obliquity: f64) -> posidonius::Particle {
    match planet_evolution {
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be non evolving to create a terrestrial/earth-like planet!"); }
    }
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let planet_radius = 1.5687291659 * posidonius::constants::R_EARTH;

    ////// Radius of gyration
    let planet_radius_of_gyration: f64 = 5.75e-01; // Earth type planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = semimajor_axis;                             // semi-major axis (in AU)
    let e: f64 = eccentricity;                               // eccentricity
    let i: f64 = 0. * posidonius::constants::DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0. * posidonius::constants::DEG2RAD;                            // argument of pericentre (degrees)
    let n: f64 = 0. * posidonius::constants::DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * posidonius::constants::DEG2RAD;                      // mean anomaly (degrees)
    p = p + n;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    let (x, y, z, vx, vy, vz) = posidonius::tools::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let planet_position = posidonius::Axes{x:x, y:y, z:z};
    let planet_velocity = posidonius::Axes{x:vx, y:vy, z:vz};

    ////// Initialization of planetary spin
    let orbital_period = posidonius::constants::TWO_PI * (((a*posidonius::constants::AU).powi(3)) / ( posidonius::constants::G_SI*( star.mass*posidonius::constants::M_SUN + planet_mass*posidonius::constants::M_SUN))).sqrt()/86400.; // in days
    let planet_rotation_period = orbital_period * (1.03); // days
    let planet_angular_frequency = posidonius::constants::TWO_PI/(planet_rotation_period); // days^-1
    //
    // Leon: Following statement seems to be a bit redundant/duplicate, see if one can reduce it.
    let planet_keplerian_orbital_elements = posidonius::tools::calculate_keplerian_orbital_elements(posidonius::constants::G*star.mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    let planet_longitude_ascending_node = 0.;
    let planet_spin = posidonius::tools::calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_longitude_ascending_node);


    /////////////////////////////////////////////////////////////////////////////////
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::DisabledModel));
    let planet_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    let planet_general_relativity = posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::Disabled);
    let planet_wind_k_factor = 0.;
    let planet_wind_rotation_saturation = 0.;
    let planet_wind = posidonius::Wind::new(posidonius::WindEffect::Disabled, planet_wind_k_factor, planet_wind_rotation_saturation);
    let planet_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut planet = posidonius::Particle::new(planet_mass, planet_radius, planet_radius_of_gyration,
                                               planet_position, planet_velocity, planet_spin);
    planet.set_tides(planet_tides);
    planet.set_rotational_flattening(planet_rotational_flattening);
    planet.set_general_relativity(planet_general_relativity);
    planet.set_wind(planet_wind);
    planet.set_disk(planet_disk);
    planet.set_evolution(planet_evolution);
    planet
}

#[allow(dead_code)]
pub fn trappist1_h_like(star: &posidonius::Particle, planet_mass: f64, planet_evolution: posidonius::EvolutionType, semimajor_axis: f64, eccentricity: f64, planet_obliquity: f64) -> posidonius::Particle {
    match planet_evolution {
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be non evolving to create a terrestrial/earth-like planet!"); }
    }
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    let planet_radius = 4810105.0000000000 / posidonius::constants::AU; // Assume same average density with Earth's

    ////// Radius of gyration
    let gyration_radius_2: f64 = 0.33916953912205672;
    let planet_radius_of_gyration: f64 = (gyration_radius_2).sqrt();

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    let a: f64 = semimajor_axis;                             // semi-major axis (in AU)
    let e: f64 = eccentricity;                               // eccentricity
    let i: f64 = 0. * posidonius::constants::DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0. * posidonius::constants::DEG2RAD;                            // argument of pericentre (degrees)
    let n: f64 = 185. * posidonius::constants::DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 290.85550 * posidonius::constants::DEG2RAD;  // mean anomaly (degrees)
    p = p + n;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    let (x, y, z, vx, vy, vz) = posidonius::tools::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    let planet_position = posidonius::Axes{x:x, y:y, z:z};
    let planet_velocity = posidonius::Axes{x:vx, y:vy, z:vz};

    ////// Initialization of planetary spin
    let orbital_period = posidonius::constants::TWO_PI * (((a*posidonius::constants::AU).powi(3)) / ( posidonius::constants::G_SI*( star.mass*posidonius::constants::M_SUN + planet_mass*posidonius::constants::M_SUN))).sqrt()/86400.; // in days
    let planet_rotation_period = orbital_period * (1.03); // days
    let planet_angular_frequency = posidonius::constants::TWO_PI/(planet_rotation_period); // days^-1
    // dbg!(orbital_period, planet_rotation_period, planet_angular_frequency);

    //
    // Leon: Following statement seems to be a bit redundant/duplicate, see if one can reduce it.
    let planet_keplerian_orbital_elements = posidonius::tools::calculate_keplerian_orbital_elements(posidonius::constants::G*star.mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    let planet_longitude_ascending_node = 0.;
    let planet_spin = posidonius::tools::calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_longitude_ascending_node);
    // dbg!(planet_spin, planet_angular_frequency, planet_inclination, planet_obliquity, planet_longitude_ascending_node);
    // panic!();

    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::DisabledModel));
    ///////////////////////////////////////////////////////////////////////////////////
    let planet_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    let planet_general_relativity = posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::Disabled);
    let planet_wind_k_factor = 0.;
    let planet_wind_rotation_saturation = 0.;
    let planet_wind = posidonius::Wind::new(posidonius::WindEffect::Disabled, planet_wind_k_factor, planet_wind_rotation_saturation);
    let planet_disk = posidonius::Disk::new(posidonius::DiskEffect::Disabled);
    let mut planet = posidonius::Particle::new(planet_mass, planet_radius, planet_radius_of_gyration,
                                               planet_position, planet_velocity, planet_spin);
    planet.set_tides(planet_tides);
    planet.set_rotational_flattening(planet_rotational_flattening);
    planet.set_general_relativity(planet_general_relativity);
    planet.set_wind(planet_wind);
    planet.set_disk(planet_disk);
    planet.set_evolution(planet_evolution);
    planet
}

#[allow(dead_code)]
pub fn earth_like(star: &posidonius::Particle, planet_mass: f64, planet_evolution: posidonius::EvolutionType, semimajor_axis: f64, eccentricity: f64, planet_obliquity: f64) -> posidonius::Particle {
    match planet_evolution {
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be non evolving to create a terrestrial/earth-like planet!"); }
    }
    ////// Planetary radius in AU (rearth in AU) Rocky planet
    // Earth-like => mass-radius relationship from Fortney 2007
    let planet_percent_rock = 0.70;
    let planet_mass_factor: f64 = planet_mass / posidonius::constants::M_EARTH;
    let planet_radius_factor : f64 = (0.0592*planet_percent_rock+0.0975) * planet_mass_factor.log10().powi(2)
        + (0.2337*planet_percent_rock+0.4938) * planet_mass_factor.log10()
        + 0.3102*planet_percent_rock+0.7932;
    //let planet_radius_factor : f64 = 0.00004307581066392812;
    let planet_radius: f64 = planet_radius_factor * posidonius::constants::R_EARTH;
    ////// Radius of gyration
    let planet_radius_of_gyration: f64 = 5.75e-01; // Earth type planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    //let a: f64 = 0.018;                             // semi-major axis (in AU)
    let a: f64 = semimajor_axis;                             // semi-major axis (in AU)
    //let e: f64 = 0.1;                               // eccentricity
    let e: f64 = eccentricity;                               // eccentricity
    let i: f64 = 5. * posidonius::constants::DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0. * posidonius::constants::DEG2RAD;                            // argument of pericentre (degrees)
    let n: f64 = 0. * posidonius::constants::DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * posidonius::constants::DEG2RAD;                      // mean anomaly (degrees)
    p = p + n;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    let (x, y, z, vx, vy, vz) = posidonius::tools::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    //////// Cartesian coordinates
    //let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    //let x: f64 = 0.1;
    //let y: f64 = 0.;
    //let z: f64 = 0.;
    //let vx: f64 = 0.;
    //let vy: f64 = core::num::Float::sqrt(gm/x); // Stable orbit
    //let vz: f64 = 0.;

    let planet_position = posidonius::Axes{x:x, y:y, z:z};
    let planet_velocity = posidonius::Axes{x:vx, y:vy, z:vz};

    ////// Initialization of planetary spin
    //let planet_obliquity: f64 = 11.459156 * posidonius::constants::DEG2RAD; // 0.2 rad
    let planet_rotation_period: f64 = 24.; // hours
    let planet_angular_frequency = posidonius::constants::TWO_PI/(planet_rotation_period/24.); // days^-1
    //
    let planet_keplerian_orbital_elements = posidonius::tools::calculate_keplerian_orbital_elements(posidonius::constants::G*star.mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    let planet_longitude_ascending_node = 0.;
    let planet_spin = posidonius::tools::calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_longitude_ascending_node);

    let planet_love_number: f64 = 0.299; // Earth
    let planet_fluid_love_number: f64 = 0.9532; // Earth
    ////// Disipation factor (sigma)
    let planet_dissipation_factor_scale: f64 = 1.;
    //// Hot Gas Giant:
    //let planet_dissipation_factor: f64 = 2.006*3.845764d4;
    //// Terrestrial:
    let k2pdelta: f64 = 2.465278e-3; // Terrestrial planets (no gas)
    let planet_dissipation_factor: f64 = 2. * posidonius::constants::K2 * k2pdelta/(3. * planet_radius.powi(5));
    /////////////////////////////////////////////////////////////////////////////////
    // Constant Time Lag
    let planet_tidal_model_params = posidonius::ConstantTimeLagParameters{dissipation_factor: planet_dissipation_factor,
        dissipation_factor_scale: planet_dissipation_factor_scale,
        love_number: planet_love_number};
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::ConstantTimeLag(planet_tidal_model_params)));
    ///////////////////////////////////////////////////////////////////////////////////
    let planet_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters{love_number: planet_fluid_love_number};
    let planet_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::OrbitingBody(posidonius::RotationalFlatteningModel::OblateSpheroid(planet_rotational_flattening_oblate_spheroid_params)));
    let planet_general_relativity = posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::OrbitingBody);
    let planet_wind_k_factor = 0.;
    let planet_wind_rotation_saturation = 0.;
    let planet_wind = posidonius::Wind::new(posidonius::WindEffect::Disabled, planet_wind_k_factor, planet_wind_rotation_saturation);
    let planet_disk = posidonius::Disk::new(posidonius::DiskEffect::OrbitingBody);
    let mut planet = posidonius::Particle::new(planet_mass, planet_radius, planet_radius_of_gyration,
                                               planet_position, planet_velocity, planet_spin);
    planet.set_tides(planet_tides);
    planet.set_rotational_flattening(planet_rotational_flattening);
    planet.set_general_relativity(planet_general_relativity);
    planet.set_wind(planet_wind);
    planet.set_disk(planet_disk);
    planet.set_evolution(planet_evolution);
    planet
}

#[allow(dead_code)]
pub fn jupiter_like(star: &posidonius::Particle, planet_mass: f64, planet_evolution: posidonius::EvolutionType, semimajor_axis: f64, eccentricity: f64, planet_obliquity: f64) -> posidonius::Particle {
    match planet_evolution {
        posidonius::EvolutionType::LeconteChabrier2013(_) => { },
        posidonius::EvolutionType::NonEvolving => { },
        _ => { panic!("Evolution type should be LeconteChabrier2013 or non evolving to create a gaseous/jupiter-like planet!"); }
    }

    // Planetary radius in AU (rearth in AU)
    let planet_radius_factor: f64 = 10.9;
    let planet_radius: f64 = planet_radius_factor * posidonius::constants::R_EARTH;
    // Radius of gyration
    let planet_radius_of_gyration: f64 = 5.04e-01; // Gas planet

    ////////// Specify initial position and velocity for a stable orbit
    ////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    //let a: f64 = 0.018;                             // semi-major axis (in AU)
    let a: f64 = semimajor_axis;                             // semi-major axis (in AU)
    //let e: f64 = 0.1;                               // eccentricity
    let e: f64 = eccentricity;                               // eccentricity
    let i: f64 = 5. * posidonius::constants::DEG2RAD;                      // inclination (degrees)
    let mut p: f64 = 0. * posidonius::constants::DEG2RAD;                            // argument of pericentre (degrees)
    let n: f64 = 0. * posidonius::constants::DEG2RAD;                      // longitude of the ascending node (degrees)
    let l: f64 = 0. * posidonius::constants::DEG2RAD;                      // mean anomaly (degrees)
    p = p + n;                          // Convert to longitude of perihelion !!
    let q = a * (1.0 - e);                          // perihelion distance
    let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    let (x, y, z, vx, vy, vz) = posidonius::tools::calculate_cartesian_coordinates(gm, q, e, i, p, n, l);

    //////// Cartesian coordinates
    //let gm: f64 = posidonius::constants::G*(planet_mass+star.mass);
    //let x: f64 = 0.1;
    //let y: f64 = 0.;
    //let z: f64 = 0.;
    //let vx: f64 = 0.;
    //let vy: f64 = core::num::Float::sqrt(gm/x); // Stable orbit
    //let vz: f64 = 0.;

    let planet_position = posidonius::Axes{x:x, y:y, z:z};
    let planet_velocity = posidonius::Axes{x:vx, y:vy, z:vz};

    ////// Initialization of planetary spin
    //let planet_obliquity: f64 = 11.459156 * posidonius::constants::DEG2RAD; // 0.2 rad
    let planet_rotation_period: f64 = 24.; // hours
    let planet_angular_frequency = posidonius::constants::TWO_PI/(planet_rotation_period/24.); // days^-1
    //
    let planet_keplerian_orbital_elements = posidonius::tools::calculate_keplerian_orbital_elements(posidonius::constants::G*star.mass*planet_mass, planet_position, planet_velocity);
    let planet_inclination = planet_keplerian_orbital_elements.3;
    let planet_longitude_ascending_node = 0.;
    let planet_spin = posidonius::tools::calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_longitude_ascending_node);

    let planet_love_number: f64 = 0.380; // Gas
    // Disipation factor (sigma)
    let planet_dissipation_factor_scale: f64 = 1.;
    // Hot Gas Giant:
    //let k2pdelta: f64 = 8.101852e-9; // Gas giant
    let k2pdelta: f64 = 2.893519e-7; // Gas giant for Jupiter: 2-3d-2 s, here in day (Leconte)
    let planet_dissipation_factor: f64 = 2. * posidonius::constants::K2 * k2pdelta/(3. * planet_radius.powi(5));
    //let planet_dissipation_factor: f64 = 2.006*3.845764d4; // Gas giant
    let planet_tidal_model_params = posidonius::ConstantTimeLagParameters{dissipation_factor: planet_dissipation_factor,
        dissipation_factor_scale: planet_dissipation_factor_scale,
        love_number: planet_love_number};
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(posidonius::TidalModel::ConstantTimeLag(planet_tidal_model_params)));
    let planet_rotational_flattening_oblate_spheroid_params = posidonius::OblateSpheroidParameters{love_number: planet_love_number};
    let planet_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::OrbitingBody(posidonius::RotationalFlatteningModel::OblateSpheroid(planet_rotational_flattening_oblate_spheroid_params)));
    let planet_general_relativity = posidonius::GeneralRelativity::new(posidonius::GeneralRelativityEffect::OrbitingBody);
    let planet_wind_k_factor = 0.;
    let planet_wind_rotation_saturation = 0.;
    let planet_wind = posidonius::Wind::new(posidonius::WindEffect::Disabled, planet_wind_k_factor, planet_wind_rotation_saturation);
    let planet_disk = posidonius::Disk::new(posidonius::DiskEffect::OrbitingBody);
    let mut planet = posidonius::Particle::new(planet_mass, planet_radius, planet_radius_of_gyration,
                                               planet_position, planet_velocity, planet_spin);
    planet.set_tides(planet_tides);
    planet.set_rotational_flattening(planet_rotational_flattening);
    planet.set_general_relativity(planet_general_relativity);
    planet.set_wind(planet_wind);
    planet.set_disk(planet_disk);
    planet.set_evolution(planet_evolution);
    planet
}
