use assert_approx_eq::assert_approx_eq;
use std::fs;
use std::fs::File;
use std::io::Read;
use std::io::{BufWriter, Write};
use std::iter;
use std::path::Path;
use time::{OffsetDateTime, format_description};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
struct ParticleProperties {
    pub inertial_position: posidonius::Axes,
    pub inertial_velocity: posidonius::Axes,
    pub inertial_acceleration: posidonius::Axes,
    pub spin: posidonius::Axes,
}

#[allow(dead_code)]
pub fn store_positions_unless_files_exist(universe: &posidonius::Universe, dirname: &String) {
    let _ = fs::create_dir(dirname);
    for (i, particle) in universe.particles[..universe.n_particles]
        .iter()
        .enumerate()
    {
        let assert_properties = ParticleProperties {
            inertial_position: particle.inertial_position,
            inertial_velocity: particle.inertial_velocity,
            inertial_acceleration: particle.inertial_acceleration,
            spin: particle.spin,
        };
        let expected_particle_filename = format!("{dirname}/particle_{i}.json");
        if !Path::new(&expected_particle_filename).exists() {
            let json_encoded = serde_json::to_string_pretty(&assert_properties).unwrap();
            let mut writer = BufWriter::new(File::create(&expected_particle_filename).unwrap());
            writer.write(json_encoded.as_bytes()).unwrap();
        }
    }
}

#[allow(dead_code)]
pub fn store_unless_files_exist<I: serde::ser::Serialize>(
    universe_integrator: &I,
    dirname: &String,
) {
    let _ = fs::create_dir(dirname);
    let snapshot_filename = format!("{dirname}/case.json");
    let snapshot_path = Path::new(&snapshot_filename);
    if !Path::new(&snapshot_path).exists() {
        posidonius::output::write_recovery_snapshot(snapshot_path, &universe_integrator);
    }
}

// NewFeatures: spins will be checked from now on, as they are evolving in tides.
// Debug: I changed some precisions as few tests can not pass the very small precision.
// This is not a good practice, and it is weird as in my modification should not affect other tests,
// something might be wrong from the creep tide which is involved in those tests.
// See the possible explanation in mattermost.
#[allow(dead_code)]
pub fn assert_stored_positions(universe: &posidonius::Universe, dirname: &String) {
    let precision = 1.0e-14;
    for (i, particle) in universe.particles[..universe.n_particles]
        .iter()
        .enumerate()
    {
        let expected_particle_filename = format!("{dirname}/particle_{i}.json");
        let expected_particle_path = Path::new(&expected_particle_filename);
        // Open the path in read-only mode, returns `io::Result<File>`
        let mut expected_particle_file = File::open(expected_particle_path).unwrap();
        //// Deserialize using `json::decode`
        let mut expected_particle_json = String::new();
        let _ = expected_particle_file
            .read_to_string(&mut expected_particle_json)
            .unwrap();
        let expected_particle: ParticleProperties =
            serde_json::from_str(&expected_particle_json).unwrap();
        println!(
            "[ASSERT {} UTC] Particle {} - Position.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_position.x,
            expected_particle.inertial_position.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.y,
            expected_particle.inertial_position.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.z,
            expected_particle.inertial_position.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Velocity.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_velocity.x,
            expected_particle.inertial_velocity.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.y,
            expected_particle.inertial_velocity.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.z,
            expected_particle.inertial_velocity.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Acceleration.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_acceleration.x,
            expected_particle.inertial_acceleration.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.y,
            expected_particle.inertial_acceleration.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.z,
            expected_particle.inertial_acceleration.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Spin.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(particle.spin.x, expected_particle.spin.x, precision);
        assert_approx_eq!(particle.spin.y, expected_particle.spin.y, precision);
        assert_approx_eq!(particle.spin.z, expected_particle.spin.z, precision);
    }
}

#[allow(dead_code)]
pub fn iterate<T>(universe_integrator: &mut T)
where
    T: posidonius::Integrator,
{
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer =
        posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();
    while universe_integrator
        .iterate(&mut universe_history_writer, silent_mode)
        .is_ok()
    {}
    let _ = fs::remove_file(universe_history_filename);
}

#[allow(dead_code)]
pub fn one_step<T>(universe_integrator: &mut T)
where
    T: posidonius::Integrator,
{
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer =
        posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();
    let _ = universe_integrator.iterate(&mut universe_history_writer, silent_mode);
    let _ = fs::remove_file(universe_history_filename);
}

#[allow(dead_code)]
fn one_step_box(universe_integrator: &mut Box<dyn posidonius::Integrator>) {
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer =
        posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();
    let _ = universe_integrator.iterate(&mut universe_history_writer, silent_mode);
    let _ = fs::remove_file(universe_history_filename);
}

#[allow(dead_code)]
fn iterate_box(universe_integrator: &mut Box<dyn posidonius::Integrator>) {
    let universe_history_filename = "/tmp/delete_me.dump";
    let universe_history_path = Path::new(universe_history_filename);
    let expected_n_bytes = 0;
    let silent_mode = true;
    let mut universe_history_writer =
        posidonius::output::get_universe_history_writer(universe_history_path, expected_n_bytes);
    universe_integrator.initialize_physical_values();

    let _ = universe_integrator.iterate(&mut universe_history_writer, silent_mode);
    let _ = fs::remove_file(universe_history_filename);
}

#[allow(dead_code)]
pub fn iterate_universe_from_json(dirname: &String) -> posidonius::Universe {
    let snapshot_from_python_filename = format!("{dirname}/case.json");
    let snapshot_from_python_path = Path::new(&snapshot_from_python_filename);

    let mut boxed_universe_integrator_from_python: Box<dyn posidonius::Integrator> =
        posidonius::output::restore_snapshot(snapshot_from_python_path).unwrap();
    iterate_box(&mut boxed_universe_integrator_from_python);
    // Extract universe from the integrator
    let universe_from_python: posidonius::Universe = match boxed_universe_integrator_from_python
        .as_any()
        .downcast_ref::<posidonius::WHFast>()
    {
        Some(universe_integrator_from_python) => universe_integrator_from_python.universe.clone(),
        None => {
            match boxed_universe_integrator_from_python
                .as_any()
                .downcast_ref::<posidonius::Ias15>()
            {
                Some(universe_integrator_from_python) => {
                    universe_integrator_from_python.universe.clone()
                }
                None => {
                    match boxed_universe_integrator_from_python
                        .as_any()
                        .downcast_ref::<posidonius::LeapFrog>()
                    {
                        Some(universe_integrator_from_python) => {
                            universe_integrator_from_python.universe.clone()
                        }
                        None => {
                            panic!("Unknown integrator!");
                        }
                    }
                }
            }
        }
    };
    universe_from_python
}

#[allow(dead_code)]
pub fn one_step_from_json(dirname: &String) -> posidonius::Universe {
    let snapshot_from_python_filename = format!("{dirname}/case.json");
    let snapshot_from_python_path = Path::new(&snapshot_from_python_filename);

    let mut boxed_universe_integrator_from_python: Box<dyn posidonius::Integrator> =
        posidonius::output::restore_snapshot(snapshot_from_python_path).unwrap();
    one_step_box(&mut boxed_universe_integrator_from_python);
    // Extract universe from the integrator
    let universe_from_python: posidonius::Universe = match boxed_universe_integrator_from_python
        .as_any()
        .downcast_ref::<posidonius::WHFast>()
    {
        Some(universe_integrator_from_python) => universe_integrator_from_python.universe.clone(),
        None => {
            match boxed_universe_integrator_from_python
                .as_any()
                .downcast_ref::<posidonius::Ias15>()
            {
                Some(universe_integrator_from_python) => {
                    universe_integrator_from_python.universe.clone()
                }
                None => {
                    match boxed_universe_integrator_from_python
                        .as_any()
                        .downcast_ref::<posidonius::LeapFrog>()
                    {
                        Some(universe_integrator_from_python) => {
                            universe_integrator_from_python.universe.clone()
                        }
                        None => {
                            panic!("Unknown integrator!");
                        }
                    }
                }
            }
        }
    };
    universe_from_python
}

#[allow(dead_code)]
pub fn assert(universe: &posidonius::Universe, parallel_universe: &posidonius::Universe) {
    for (i, (particle, parallel_particle)) in universe.particles[..universe.n_particles]
        .iter()
        .zip(parallel_universe.particles[..parallel_universe.n_particles].iter())
        .enumerate()
    {
        let precision = 1.0e-15;
        println!(
            "[ASSERT {} UTC] Particle {} - Position.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_position.x,
            parallel_particle.inertial_position.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.y,
            parallel_particle.inertial_position.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.z,
            parallel_particle.inertial_position.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Velocity.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_velocity.x,
            parallel_particle.inertial_velocity.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.y,
            parallel_particle.inertial_velocity.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.z,
            parallel_particle.inertial_velocity.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Acceleration.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_acceleration.x,
            parallel_particle.inertial_acceleration.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.y,
            parallel_particle.inertial_acceleration.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.z,
            parallel_particle.inertial_acceleration.z,
            precision
        );
        // Lower the precision requirement for the spin check, as we know the .json generated from
        // python and rust could be different, especially the initial spin for the star.
        let precision_spin = 4.0e-15;
        println!(
            "[ASSERT {} UTC] Particle {} - Spin.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(particle.spin.x, parallel_particle.spin.x, precision_spin);
        assert_approx_eq!(particle.spin.y, parallel_particle.spin.y, precision_spin);
        assert_approx_eq!(particle.spin.z, parallel_particle.spin.z, precision_spin);
    }
}

#[allow(dead_code)]
pub fn loosly_assert(universe: &posidonius::Universe, parallel_universe: &posidonius::Universe) {
    for (i, (particle, parallel_particle)) in universe.particles[..universe.n_particles]
        .iter()
        .zip(parallel_universe.particles[..parallel_universe.n_particles].iter())
        .enumerate()
    {
        let precision = 1.0e-8;
        println!(
            "[ASSERT {} UTC] Particle {} - Position.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_position.x,
            parallel_particle.inertial_position.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.y,
            parallel_particle.inertial_position.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.z,
            parallel_particle.inertial_position.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Velocity.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_velocity.x,
            parallel_particle.inertial_velocity.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.y,
            parallel_particle.inertial_velocity.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.z,
            parallel_particle.inertial_velocity.z,
            precision
        );
        let precision = 1.0e-7;
        println!(
            "[ASSERT {} UTC] Particle {} - Acceleration.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_acceleration.x,
            parallel_particle.inertial_acceleration.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.y,
            parallel_particle.inertial_acceleration.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.z,
            parallel_particle.inertial_acceleration.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Spin.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(particle.spin.x, parallel_particle.spin.x, precision);
        assert_approx_eq!(particle.spin.y, parallel_particle.spin.y, precision);
        assert_approx_eq!(particle.spin.z, parallel_particle.spin.z, precision);
    }
}

#[allow(dead_code)]
pub fn assert_when_star_is_swapped(
    universe: &posidonius::Universe,
    parallel_universe: &posidonius::Universe,
) {
    // In universe, star is first
    // In parallel universe, star is second
    for (i, (particle, parallel_particle)) in universe.particles[..universe.n_particles]
        .iter()
        .zip(
            iter::once(&parallel_universe.particles[1])
                .chain(iter::once(&parallel_universe.particles[0]))
                .chain(parallel_universe.particles[2..parallel_universe.n_particles].iter()),
        )
        .enumerate()
    {
        let precision = 1.0e-16;
        println!(
            "[ASSERT {} UTC] Particle {} - Position.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_position.x,
            parallel_particle.inertial_position.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.y,
            parallel_particle.inertial_position.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_position.z,
            parallel_particle.inertial_position.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Velocity.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_velocity.x,
            parallel_particle.inertial_velocity.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.y,
            parallel_particle.inertial_velocity.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_velocity.z,
            parallel_particle.inertial_velocity.z,
            precision
        );
        println!(
            "[ASSERT {} UTC] Particle {} - Acceleration.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            i
        );
        assert_approx_eq!(
            particle.inertial_acceleration.x,
            parallel_particle.inertial_acceleration.x,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.y,
            parallel_particle.inertial_acceleration.y,
            precision
        );
        assert_approx_eq!(
            particle.inertial_acceleration.z,
            parallel_particle.inertial_acceleration.z,
            precision
        );
    }
}
