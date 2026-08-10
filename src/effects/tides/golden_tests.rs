// Golden regression tests for the tidal-force pipeline.
//
// These tests freeze the numerical output of the tidal pipeline for a set of
// representative configurations, so that structural refactors of the tides
// code can be verified to preserve the physics. The golden values were
// generated on commit 491c643 ("cleanup, refactor, optimise (#3)") by running
// these same scenarios and capturing the outputs.
//
// The scenarios are self-contained: no external input files (evolution tables
// or love-number spectra) are needed. Kaula spectra are synthetic but smooth
// and physically shaped; the goal is regression, not realism.
//
// The pipeline call sequence mirrors Universe::initialize +
// Universe::calculate_additional_effects for a tides-only simulation where
// the tidal host is the most massive body.

use std::collections::HashMap;

use super::constant_time_lag;
use super::kaula::{KaulaParameters, LoveNumber};
use super::{TidalModel, Tides, TidesEffect};
use crate::constants::{DAY, K2, M_EARTH, R_SUN, TWO_PI};
use crate::{Axes, EvolutionType, Particle};

// Relative tolerance for golden comparisons: loose enough to allow
// floating-point re-association from refactors, tight enough to catch any
// physics change.
const REL_TOLERANCE: f64 = 1.0e-12;
// Values below this magnitude are compared absolutely (they are zero-ish).
const ABS_FLOOR: f64 = 1.0e-30;

fn assert_all_close(observed: &[f64], golden: &[f64], config: &str) {
    assert_eq!(
        observed.len(),
        golden.len(),
        "[{config}] output length changed: {} vs golden {}",
        observed.len(),
        golden.len()
    );
    for (i, (&o, &g)) in observed.iter().zip(golden.iter()).enumerate() {
        let scale = abs!(g);
        if scale < ABS_FLOOR {
            assert!(
                abs!(o) < ABS_FLOOR,
                "[{config}] element {i}: expected ~0 ({g:e}), got {o:e}"
            );
        } else {
            let relative_difference = abs!(o - g) / scale;
            assert!(
                relative_difference < REL_TOLERANCE,
                "[{config}] element {i}: golden {g:.17e}, observed {o:.17e}, rel diff {relative_difference:e}"
            );
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Particle builders (self-contained, tides-only; all other effects disabled)
////////////////////////////////////////////////////////////////////////////////

fn star_spin_aligned() -> Axes {
    let rotation_period_days = 8.0 / 24.0; // 8 hours
    Axes::from(0., 0., TWO_PI / rotation_period_days)
}

fn star_spin_tilted() -> Axes {
    let rotation_period_days = 8.0 / 24.0;
    let spin_z = TWO_PI / rotation_period_days;
    // Significant x component to force the 3D (oblique) Kaula branch
    Axes::from(0.3 * spin_z, 0., spin_z)
}

fn ctl_star_params() -> constant_time_lag::ConstantTimeLagParameters {
    constant_time_lag::ConstantTimeLagParameters {
        dissipation_factor: 4.992 * 3.845764e-2,
        dissipation_factor_scale: 1.0,
        love_number: 0.03,
    }
}

fn ctl_planet_params() -> constant_time_lag::ConstantTimeLagParameters {
    constant_time_lag::ConstantTimeLagParameters {
        dissipation_factor: 2.006 * 3.845764e4,
        dissipation_factor_scale: 1.0,
        love_number: 0.305,
    }
}

fn make_star(effect: TidesEffect, spin: Axes, evolution: EvolutionType) -> Particle {
    let mass = 1.0; // M_SUN
    let radius = R_SUN;
    let radius_of_gyration = 2.43e-1;
    let mut star = Particle::new(
        mass,
        radius,
        radius_of_gyration,
        Axes::new(),
        Axes::new(),
        spin,
    );
    star.id = 0;
    star.tides = Tides::new(effect);
    star.evolution = evolution;
    // Host at origin: heliocentric distance and radial velocity are zero.
    star
}

fn make_planet(effect: TidesEffect, id: usize, semi_major_axis: f64) -> Particle {
    let mass = M_EARTH;
    let radius = 1.0 * 4.26352e-5; // ~ Earth radius in AU
    let radius_of_gyration = 5.75e-1;
    let position = Axes::from(semi_major_axis, 0., 0.);
    // Slightly sub-circular velocity for a mildly eccentric orbit (e ~ 0.04)
    let circular_velocity = sqrt!(K2 * (1.0 + mass) / semi_major_axis);
    let velocity = Axes::from(0., 0.98 * circular_velocity, 0.);
    // Planet spin: 24h period, aligned with z
    let spin = Axes::from(0., 0., TWO_PI / 1.0);
    let mut planet = Particle::new(mass, radius, radius_of_gyration, position, velocity, spin);
    planet.id = id;
    planet.tides = Tides::new(effect);
    planet.evolution = EvolutionType::NonEvolving;
    planet.heliocentric_distance = sqrt!(
        position.x * position.x + position.y * position.y + position.z * position.z
    );
    planet.heliocentric_radial_velocity = (position.x * velocity.x
        + position.y * velocity.y
        + position.z * velocity.z)
        / planet.heliocentric_distance;
    planet
}

// Synthetic but smooth love-number spectrum: sorted excitation frequencies
// with an odd sigmoid-shaped imaginary part and a slowly varying real part.
fn synthetic_love_number(stellar_spectrum_spin_rate: Option<f64>) -> LoveNumber {
    let mut frequency = [0.; 1024];
    let mut real_part = [0.; 1024];
    let mut imaginary_part = [0.; 1024];
    let width = 2.0e-4; // rad/s, covers typical wk2 magnitudes for close-in orbits
    for i in 0..1024 {
        let f = -width + 2.0 * width * (i as f64) / 1023.0;
        frequency[i] = f;
        imaginary_part[i] = 1.0e-2 * f / (1.0e-4 + abs!(f));
        real_part[i] = -0.3 + 0.01 * (f * 1.0e4).tanh();
    }
    LoveNumber::new_from(
        frequency,
        real_part,
        imaginary_part,
        stellar_spectrum_spin_rate,
    )
}

fn kaula_params(stellar_spectrum_spin_rate: Option<f64>) -> KaulaParameters {
    KaulaParameters {
        love_numbers: synthetic_love_number(stellar_spectrum_spin_rate),
        tidal_force: Axes::new(),
        polynomials: super::Polynomials::default(),
    }
}

////////////////////////////////////////////////////////////////////////////////
// Pipeline runner: mirrors Universe::initialize + calculate_additional_effects
// (tides-only, host == most massive) followed by calculate_denergy_dt.
////////////////////////////////////////////////////////////////////////////////

fn run_tidal_pipeline(star: &mut Particle, planets: &mut [Particle]) -> Vec<f64> {
    let more: &mut [Particle] = &mut [];
    let mut pair_dependent_scaled_dissipation_factor: HashMap<usize, f64> = HashMap::new();

    super::copy_heliocentric_coordinates(star, planets, more);
    super::initialize(star, planets, more);
    constant_time_lag::calculate_pair_dependent_scaled_dissipation_factors(
        star,
        planets,
        more,
        &mut pair_dependent_scaled_dissipation_factor,
    );
    constant_time_lag::calculate_orthogonal_component_of_the_tidal_force(
        star,
        planets,
        more,
        &mut pair_dependent_scaled_dissipation_factor,
    );
    super::calculate_creep_coplanar_shapes(star, planets, more);
    constant_time_lag::calculate_radial_component_of_the_tidal_force(
        star,
        planets,
        more,
        &mut pair_dependent_scaled_dissipation_factor,
    );
    super::calculate_tidal_acceleration(star, planets, more);
    super::calculate_dangular_momentum_dt_due_to_tides(star, planets, more);
    super::calculate_denergy_dt(planets, more);

    // Flatten all physically meaningful outputs into a single vector.
    let mut output = Vec::new();
    output.push(star.tides.parameters.output.acceleration.x);
    output.push(star.tides.parameters.output.acceleration.y);
    output.push(star.tides.parameters.output.acceleration.z);
    output.push(star.tides.parameters.output.dangular_momentum_dt.x);
    output.push(star.tides.parameters.output.dangular_momentum_dt.y);
    output.push(star.tides.parameters.output.dangular_momentum_dt.z);
    for planet in planets.iter() {
        output.push(planet.tides.parameters.output.acceleration.x);
        output.push(planet.tides.parameters.output.acceleration.y);
        output.push(planet.tides.parameters.output.acceleration.z);
        output.push(planet.tides.parameters.output.dangular_momentum_dt.x);
        output.push(planet.tides.parameters.output.dangular_momentum_dt.y);
        output.push(planet.tides.parameters.output.dangular_momentum_dt.z);
        output.push(planet.tides.parameters.internal.denergy_dt);
    }
    // Record the pair-dependent sigma overrides (sorted by key for stability).
    let mut keys: Vec<usize> = pair_dependent_scaled_dissipation_factor
        .keys()
        .copied()
        .collect();
    keys.sort_unstable();
    for key in keys {
        output.push(key as f64);
        output.push(pair_dependent_scaled_dissipation_factor[&key]);
    }
    output
}

#[allow(dead_code)]
fn print_golden(config: &str, output: &[f64]) {
    // Print in a paste-ready format for (re)capturing golden values.
    println!("const GOLDEN_{}: &[f64] = &[", config.to_uppercase());
    for value in output {
        println!("    {value:.17e},");
    }
    println!("];");
}

////////////////////////////////////////////////////////////////////////////////
// Golden values, generated on commit 491c643 (pre-restructure physics).
//
// Layout per configuration:
//   [0..3]  star acceleration (x, y, z)
//   [3..6]  star dangular_momentum_dt (x, y, z)
//   then per planet: acceleration (x, y, z), dangular_momentum_dt (x, y, z),
//   denergy_dt — 7 values each
//   then, if any, pair-dependent sigma map entries as (key, value) pairs
////////////////////////////////////////////////////////////////////////////////

const GOLDEN_CTL_STAR_CTL_PLANET_EQUILIBRIUM: &[f64] = &[
    4.93284483305164412e-17,
    -8.18606782694432665e-25,
    -0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    -4.09282323372251416e-26,
    -1.64237120684235051e-11,
    2.72551895533972755e-19,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.10679749649641472e-30,
    1.00659900148855304e-29,
];

const GOLDEN_CTL_STAR_CTL_PLANET_DYNAMICAL: &[f64] = &[
    4.93284483305164412e-17,
    -1.20088429072989347e-23,
    -0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    -6.00440038567450292e-25,
    -1.64237120684235051e-11,
    3.99829682180361044e-18,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.10679749649641472e-30,
    1.00659900148855304e-29,
    1.00000000000000000e0,
    2.81646178167494199e0,
];

const GOLDEN_CTL_STAR_TWO_CTL_PLANETS: &[f64] = &[
    5.01341790441843755e-17,
    -8.32657943003470494e-25,
    -0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    -4.21927633881833587e-26,
    -1.64237120684235051e-11,
    2.72551895533972755e-19,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.10679749649641472e-30,
    1.00659900148855304e-29,
    -2.68264859200548384e-13,
    4.67827833538672429e-21,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -7.33768551907041022e-32,
    4.15302100825151486e-31,
];

const GOLDEN_DISABLED_MODEL_STAR_CTL_PLANET: &[f64] = &[
    4.86617901775889490e-17,
    -4.21359499299282929e-29,
    -0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    -1.62017508691093397e-11,
    1.40289981298790782e-23,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.10679749649641472e-30,
    1.00659900148855304e-29,
];

const GOLDEN_KAULA_STAR_CTL_PLANET_OBLIQUE: &[f64] = &[
    6.50334400940568501e-17,
    2.63575606909435520e-20,
    3.51734812611644071e-21,
    0.00000000000000000e0,
    -3.21956518137532754e-23,
    -2.27728466557103443e-21,
    -2.16526270554327836e-11,
    -8.77564574327494313e-15,
    -1.17108716821339727e-15,
    -0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.10679749649641472e-30,
    1.00659900148855304e-29,
];

const GOLDEN_KAULA_STAR_KAULA_PLANET_ALIGNED: &[f64] = &[
    5.38447225712696049e-17,
    -9.53184992316114282e-19,
    0.00000000000000000e0,
    0.00000000000000000e0,
    0.00000000000000000e0,
    -3.39656772689538733e-20,
    -1.79273877416410838e-11,
    3.17359178964031149e-13,
    -0.00000000000000000e0,
    0.00000000000000000e0,
    -0.00000000000000000e0,
    -2.63788241808948646e-20,
    1.46652478084938967e-20,
];

////////////////////////////////////////////////////////////////////////////////
// Scenario tests
////////////////////////////////////////////////////////////////////////////////

#[test]
fn golden_ctl_star_ctl_planet_equilibrium() {
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(ctl_star_params())),
        star_spin_aligned(),
        EvolutionType::NonEvolving,
    );
    let mut planets = [make_planet(
        TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
        1,
        0.05,
    )];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_CTL_STAR_CTL_PLANET_EQUILIBRIUM, "ctl_star_ctl_planet_equilibrium");
}

#[test]
fn golden_ctl_star_ctl_planet_dynamical() {
    // Star with an evolution type that (in the pre-refactor code) activates
    // the frequency-averaged dynamical tide via the pair-dependent sigma map.
    // The lag angle is set manually to a representative value; the Evolver
    // itself (which needs external files) is not involved.
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(ctl_star_params())),
        star_spin_aligned(),
        EvolutionType::BolmontMathis2016(1.0),
    );
    star.tides.parameters.internal.lag_angle = 5.0e-7;
    let mut planets = [make_planet(
        TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
        1,
        0.05,
    )];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_CTL_STAR_CTL_PLANET_DYNAMICAL, "ctl_star_ctl_planet_dynamical");
}

#[test]
fn golden_kaula_star_ctl_planet_oblique() {
    // Mixed configuration: Kaula stellar tide (3D branch via tilted spin),
    // CTL planetary tide.
    let spin = star_spin_tilted();
    let stellar_spectrum_spin_rate = sqrt!(
        spin.x * spin.x + spin.y * spin.y + spin.z * spin.z
    ) / DAY; // rad/s, matching the runtime unit used by kaula.rs
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::Kaula(kaula_params(Some(
            stellar_spectrum_spin_rate,
        )))),
        spin,
        EvolutionType::NonEvolving,
    );
    let mut planets = [make_planet(
        TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
        1,
        0.05,
    )];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_KAULA_STAR_CTL_PLANET_OBLIQUE, "kaula_star_ctl_planet_oblique");
}

#[test]
fn golden_kaula_star_kaula_planet_aligned() {
    // Pure Kaula configuration, 2D branch (aligned spins).
    let spin = star_spin_aligned();
    let stellar_spectrum_spin_rate = spin.z / DAY;
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::Kaula(kaula_params(Some(
            stellar_spectrum_spin_rate,
        )))),
        spin,
        EvolutionType::NonEvolving,
    );
    let mut planets = [make_planet(
        TidesEffect::OrbitingBody(TidalModel::Kaula(kaula_params(None))),
        1,
        0.05,
    )];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_KAULA_STAR_KAULA_PLANET_ALIGNED, "kaula_star_kaula_planet_aligned");
}

#[test]
fn golden_ctl_star_two_ctl_planets() {
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::ConstantTimeLag(ctl_star_params())),
        star_spin_aligned(),
        EvolutionType::NonEvolving,
    );
    let mut planets = [
        make_planet(
            TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
            1,
            0.05,
        ),
        make_planet(
            TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
            2,
            0.09,
        ),
    ];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_CTL_STAR_TWO_CTL_PLANETS, "ctl_star_two_ctl_planets");
}

#[test]
fn golden_disabled_model_star_ctl_planet() {
    // Planetary tide only: the star participates as tidal host but carries no
    // tidal model of its own (the documented DisabledModel workaround).
    let mut star = make_star(
        TidesEffect::CentralBody(TidalModel::DisabledModel),
        star_spin_aligned(),
        EvolutionType::NonEvolving,
    );
    let mut planets = [make_planet(
        TidesEffect::OrbitingBody(TidalModel::ConstantTimeLag(ctl_planet_params())),
        1,
        0.05,
    )];
    let output = run_tidal_pipeline(&mut star, &mut planets);
    assert_all_close(&output, GOLDEN_DISABLED_MODEL_STAR_CTL_PLANET, "disabled_model_star_ctl_planet");
}
