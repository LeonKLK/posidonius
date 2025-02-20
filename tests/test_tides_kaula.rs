#[macro_use]
extern crate serde_derive;

mod common;
use posidonius::constants::{DAY, TWO_PI};
use std::path::Path;

////////////////////////////////////////////////////////////////////////////////
fn enabled_star_tides_case() -> posidonius::WHFast {
    // simulation_properties provides a basic time-related setup
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) =
        common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Disabled; // Disable GR

    let star_mass: f64 = 1.; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    //let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution = posidonius::EvolutionType::NonEvolving;

    // For `star_evolution = posidonius::EvolutionType::Leconte2011(star_mass)`,
    // it provides the initial spin for the star
    let star = common::stars::solar_like_for_kaula_1(
        star_mass,
        star_evolution,
        general_relativity_implementation,
    );
    let mut particles = vec![star];
    particles.extend(common::planets::basic_planetary_configuration_for_stellar_kaula(&star));

    //////////////////////////////////////////////////////////////////////////////////
    // Leon: We are able to simply switch off the tides now, I saved this section of the code for
    // the reference of future developers, as I do not know much on the creep model.
    // `check_uniform_viscosity_coefficient` seems to be related to creep tide
    // src/particles/particle.rs, and creep model does not concern us now (6/11).
    // Disable all first or `check_uniform_viscosity_coefficient` will fail
    // let disabled_tides = posidonius::Tides::new(posidonius::TidesEffect::Disabled);
    // // let disabled_rotational_flattening = posidonius::RotationalFlattening::new(posidonius::RotationalFlatteningEffect::Disabled);
    // for particle in particles.iter_mut() {
    //     particle.set_tides(disabled_tides);
    //     // particle.set_rotational_flattening(disabled_rotational_flattening);
    // }
    //////////////////////////////////////////////////////////////////////////////////
    // Leon: Tides are set to "True" but we set Tides::disabled for all bodies' in solar_like/super_earth_like.
    // In this way the initialization functions can be more general.
    // Below we set back the tides according to the test case.
    // Load star data
    let star_file = "./input/love_numbers/Aurelie_Stellar/alpha0.51600_P1p2_Ek1p5em6.txt";
    // If we are simulating stellar tide, then we need the value of the spectrum spin rate,
    // which was an input value of stellar spin rate from the computation of the spectrum.
    let stellar_tide = true;
    let spectrum_stellar_spin_period = 1.2 * DAY;
    let spectrum_spin_rate = TWO_PI / spectrum_stellar_spin_period;
    let star_tidal_model_params =
        common::load_kaula_parameters(star_file, stellar_tide, spectrum_spin_rate).unwrap();
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::Kaula(star_tidal_model_params),
    ));

    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(
        posidonius::TidalModel::DisabledModel,
    ));

    if let Some((star, planets)) = particles.split_first_mut() {
        star.set_tides(star_tides);
        for planet in planets.iter_mut() {
            planet.set_tides(planet_tides);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    // let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    // let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    // let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    // let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);

    posidonius::WHFast::new(
        time_step,
        recovery_snapshot_period,
        historic_snapshot_period,
        universe,
        alternative_coordinates_type,
    )
}

#[test]
fn enabled_star_tides_rust() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_star_tides"
    );
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = enabled_star_tides_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(
        &universe_integrator.universe,
        &rust_data_dirname,
    );
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn star_tides_rust_vs_python() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_star_tides"
    );
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);

    let mut universe_integrator = enabled_star_tides_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::one_step(&mut universe_integrator);
    let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
    common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
    let parallel_universe_from_json =
        common::universe::iterate_universe_from_json(&python_data_dirname);
    let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
    common::universe::assert(&universe_from_json, &parallel_universe_from_json);
}

////////////////////////////////////////////////////////////////////////////////
fn enabled_planet_tides_case() -> posidonius::WHFast {
    // simulation_properties provides a basic time-related setup
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) =
        common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    // let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Disabled; // Disable GR

    let star_mass: f64 = 0.0898; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    // let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution = posidonius::EvolutionType::NonEvolving;

    // For `star_evolution = posidonius::EvolutionType::Leconte2011(star_mass)`,
    // it provides the initial spin for the star
    let star = common::stars::red_dwarf_like_for_kaula(
        star_mass,
        star_evolution,
        general_relativity_implementation,
    );

    let mut particles = vec![star];
    particles.extend(common::planets::basic_planetary_configuration_for_planetary_kaula(&star));

    //////////////////////////////////////////////////////////////////////////////////
    // let star_tidal_model_params = posidonius::ConstantTimeLagParameters{dissipation_factor: 0.,
    //     dissipation_factor_scale: 0.,
    //     love_number: 0.};
    // let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(posidonius::TidalModel::ConstantTimeLag(star_tidal_model_params)));
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::DisabledModel,
    ));

    // Load planet data
    let planet_file = "./input/love_numbers/Results_Trappist1_h_Fe_90_Si_2_170K_visco_freq_Imk2_posidonius_for_testing.txt";
    // If we are simulating stellar tide, then we need the value of the spectrum spin rate,
    // which was an input value of stellar spin rate from the computation of the spectrum.
    let stellar_tide = false;
    let spectrum_spin_rate = 0.0;
    let planet_tidal_model_params =
        common::load_kaula_parameters(planet_file, stellar_tide, spectrum_spin_rate).unwrap();
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(
        posidonius::TidalModel::Kaula(planet_tidal_model_params),
    ));
    //
    if let Some((star, planets)) = particles.split_first_mut() {
        star.set_tides(star_tides);
        for planet in planets.iter_mut() {
            planet.set_tides(planet_tides);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    // let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    // let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);

    posidonius::WHFast::new(
        time_step,
        recovery_snapshot_period,
        historic_snapshot_period,
        universe,
        alternative_coordinates_type,
    )
}

#[test]
fn enabled_planet_tides_rust() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_planet_tides"
    );
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);

    let mut universe_integrator = enabled_planet_tides_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(
        &universe_integrator.universe,
        &rust_data_dirname,
    );
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}
//
#[test]
fn planet_tides_rust_vs_python() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_planet_tides"
    );
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = enabled_planet_tides_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::one_step(&mut universe_integrator);
    let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
    common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
    let parallel_universe_from_json =
        common::universe::iterate_universe_from_json(&python_data_dirname);
    let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
    common::universe::assert(&universe_from_json, &parallel_universe_from_json);
}

////////////////////////////////////////////////////////////////////////////////

fn enabled_both_tides_case() -> posidonius::WHFast {
    let (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period) =
        common::simulation_properties();
    let consider_effects = posidonius::ConsiderEffects {
        tides: true,
        rotational_flattening: false,
        general_relativity: false,
        disk: false,
        wind: false,
        evolution: false,
    };
    // let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Kidder1995; // Mercury-T
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Anderson1975; // REBOUNDx GR
    //let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Newhall1983; // REBOUNDx GR full
    let general_relativity_implementation = posidonius::GeneralRelativityImplementation::Disabled; // Disable GR

    let star_mass: f64 = 1.; // Solar masses
    //let star_evolution = posidonius::EvolutionType::Baraffe1998(star_mass); // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    //let star_evolution = posidonius::EvolutionType::Baraffe2015(star_mass); // mass = 0.01 .. 1.40
    // let star_evolution = posidonius::EvolutionType::Leconte2011(star_mass); // BrownDwarf (mass = 0.01 .. 0.08)
    //let star_evolution = posidonius::EvolutionType::BolmontMathis2016(star_mass); // SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::GalletBolmont2017(star_mass); // SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
    //let star_evolution = posidonius::EvolutionType::LeconteChabrier2013; // Jupiter
    let star_evolution = posidonius::EvolutionType::NonEvolving;
    let star = common::stars::solar_like_for_kaula_2(
        star_mass,
        star_evolution,
        general_relativity_implementation,
    );

    let mut particles = vec![star];
    particles.extend(common::planets::basic_planetary_configuration_for_planetary_kaula(&star));
    //////////////////////////////////////////////////////////////////////////////////
    // Load star data
    let star_file = "./input/love_numbers/Aurelie_Stellar/alpha0.51600_P1p2_Ek1p5em6.txt";
    // If we are simulating stellar tide, then we need the value of the spectrum spin rate,
    // which was an input value of stellar spin rate from the computation of the spectrum.
    let stellar_tide = true;
    let spectrum_stellar_spin_period = 1.2 * DAY;
    let spectrum_spin_rate = TWO_PI / spectrum_stellar_spin_period;
    let star_tidal_model_params =
        common::load_kaula_parameters(star_file, stellar_tide, spectrum_spin_rate).unwrap();
    let star_tides = posidonius::Tides::new(posidonius::TidesEffect::CentralBody(
        posidonius::TidalModel::Kaula(star_tidal_model_params),
    ));
    // Load planet data
    let planet_file = "./input/love_numbers/Results_Trappist1_h_Fe_90_Si_2_170K_visco_freq_Imk2_posidonius_for_testing.txt";
    let stellar_tide = false;
    let spectrum_spin_rate = 0.0;
    let planet_tidal_model_params =
        common::load_kaula_parameters(planet_file, stellar_tide, spectrum_spin_rate).unwrap();
    let planet_tides = posidonius::Tides::new(posidonius::TidesEffect::OrbitingBody(
        posidonius::TidalModel::Kaula(planet_tidal_model_params),
    ));

    dbg!(
        planet_tidal_model_params.stellar_tide,
        planet_tidal_model_params.spectrum_spin_rate
    );

    //
    if let Some((star, planets)) = particles.split_first_mut() {
        star.set_tides(star_tides);
        for planet in planets.iter_mut() {
            planet.set_tides(planet_tides);
        }
    }
    //////////////////////////////////////////////////////////////////////////////////
    let universe = posidonius::Universe::new(initial_time, time_limit, particles, consider_effects);

    // let alternative_coordinates_type = posidonius::whfast::CoordinatesType::Jacobi;
    let alternative_coordinates_type = posidonius::whfast::CoordinatesType::DemocraticHeliocentric;
    //let alternative_coordinates_type = posidonius::whfast::CoordinatesType::WHDS;
    //let universe_integrator = posidonius::LeapFrog::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);
    // let universe_integrator = posidonius::Ias15::new(time_step, recovery_snapshot_period, historic_snapshot_period, universe);

    posidonius::WHFast::new(
        time_step,
        recovery_snapshot_period,
        historic_snapshot_period,
        universe,
        alternative_coordinates_type,
    )
}
//
#[test]
fn enabled_both_tides_rust() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_both_tides"
    );
    let (rust_data_dirname, _python_data_dirname) = common::get_data_dirname(&test_name);
    //
    let mut universe_integrator = enabled_both_tides_case();
    common::universe::iterate(&mut universe_integrator);
    common::universe::store_positions_unless_files_exist(
        &universe_integrator.universe,
        &rust_data_dirname,
    );
    common::universe::assert_stored_positions(&universe_integrator.universe, &rust_data_dirname);
}

#[test]
fn enabled_tides_rust_vs_python() {
    let test_name = format!(
        "{}-{}",
        Path::new(file!()).file_stem().unwrap().to_str().unwrap(),
        "enabled_both_tides"
    );
    let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
    let mut universe_integrator = enabled_both_tides_case();
    common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
    common::universe::one_step(&mut universe_integrator);
    let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
    common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
    let parallel_universe_from_json =
        common::universe::iterate_universe_from_json(&python_data_dirname);
    let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
    common::universe::assert(&universe_from_json, &parallel_universe_from_json);
}

// TODO: the json generated by python and rust is different at the array of `roche_radiuses`
// , this doesn't seem to affect us for the moment but for integrity we should take a took.
// #[test]
// fn enabled_tides_rust_vs_python_roche_radiuses_debug() {
//     let test_name = format!("{}-{}", Path::new(file!()).file_stem().unwrap().to_str().unwrap(), "enabled_both_tides");
//     let (rust_data_dirname, python_data_dirname) = common::get_data_dirname(&test_name);
//     let mut universe_integrator = enabled_both_tides_case();
//     common::universe::store_unless_files_exist(&universe_integrator, &rust_data_dirname); // Store just in case we want to inspect it/compare it to the python generated JSON
//
//     let snapshot_from_python_filename = format!("{0}/case.json", &python_data_dirname);
//     let snapshot_from_python_path = Path::new(&snapshot_from_python_filename);
//
//     let mut python_integrator_dyn = posidonius::output::restore_snapshot(&snapshot_from_python_path).unwrap();
//     let python_integrator= python_integrator_dyn.as_any().downcast_ref::<posidonius::WHFast>().unwrap();
//     assert_eq!(&universe_integrator.universe, &python_integrator.universe);
//
//     common::universe::one_step(&mut universe_integrator);
//     let parallel_universe_from_json = common::universe::one_step_from_json(&python_data_dirname);
//     common::universe::loosly_assert(&universe_integrator.universe, &parallel_universe_from_json); // Built in situ vs Restored from JSON: it requires larger precision and only one step
//     let parallel_universe_from_json = common::universe::iterate_universe_from_json(&python_data_dirname);
//     let universe_from_json = common::universe::iterate_universe_from_json(&rust_data_dirname);
//     common::universe::assert(&universe_from_json, &parallel_universe_from_json);
// }

////////////////////////////////////////////////////////////////////////////////
