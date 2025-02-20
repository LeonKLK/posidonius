use super::super::Integrator;
use super::super::constants::MIN_ORBITAL_PERIOD_TIME_STEP_RATIO;
use super::super::particles::Reference;
use super::super::particles::Universe;
use super::super::tools::calculate_keplerian_orbital_elements;
use super::super::{Axes, TidalModel, TidesEffect};
use bincode;
use serde::Serialize;
use serde_json;
use std::fs;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::{BufReader, Read};
use std::io::{BufWriter, Write};
use std::path::Path;
use time;
use time::{OffsetDateTime, format_description};

pub use super::ias15::*;
pub use super::leapfrog::*;
pub use super::whfast::*;

////////////////////////////////////////////////////////////////////////////////
//- Dump and restore functions
////////////////////////////////////////////////////////////////////////////////

pub fn write_recovery_snapshot<I: Serialize>(snapshot_path: &Path, universe_integrator: &I) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    if snapshot_path.exists() {
        //// Backup (keep only 1 per hour thanks to filename collision)
        //let new_extension = format!("{0}.bin", OffsetDateTime::now_utc().format(&format_description::parse("[year][month][day]T[hour]h)).unwrap());
        // Backup (keep only 1 every 12 hours thanks to filename collision)
        let new_extension = format!(
            "{0}.bin",
            OffsetDateTime::now_utc()
                .format(&format_description::parse("[year][month][day]T[period]").unwrap())
                .unwrap()
        );
        fs::rename(snapshot_path, snapshot_path.with_extension(new_extension)).unwrap();
    }

    // 1.- Serialize the integrator to be able to resume if the simulation is interrupted
    let mut writer = BufWriter::new(File::create(snapshot_path).unwrap());

    if snapshot_path.extension().unwrap() == "json" {
        let json_encoded = serde_json::to_string_pretty(&universe_integrator).unwrap();
        let _ = writer.write(json_encoded.as_bytes()).unwrap();
    } else {
        // Binary
        bincode::serialize_into(&mut writer, &universe_integrator).unwrap(); // bin
    }
}

pub fn n_bytes_per_particle_in_historic_snapshot() -> u64 {
    let n_stored_fields: u64 = 20;

    8 + 8 + 4 + 8 * (n_stored_fields - 3)
}

pub fn get_universe_history_writer(
    universe_history_path: &Path,
    expected_n_bytes: u64,
) -> BufWriter<File> {
    // If history snapshot exists:
    // - It validates that it contains all the expected data
    // - If it contains more, the extra data is erased

    // We create file options to write
    let mut options_bin = OpenOptions::new();
    options_bin.create(true).write(true).append(true);

    let universe_history_file = match options_bin.open(universe_history_path) {
        Ok(f) => f,
        Err(e) => panic!(
            "[PANIC {} UTC] File error: {}",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap(),
            e
        ),
    };

    let metadata = universe_history_file.metadata().unwrap();
    let current_n_bytes = metadata.len();
    assert!(
        (current_n_bytes >= expected_n_bytes),
        "[PANIC {} UTC] Historic snapshots do not contain all the expected history ({} bytes) as indicated by the recovery snapshot ({} bytes)",
        OffsetDateTime::now_utc()
            .format(
                &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                    .unwrap()
            )
            .unwrap(),
        current_n_bytes,
        expected_n_bytes
    );

    // Keep only historic data that saved until the restored snapshot (if there it is the case)
    universe_history_file.set_len(expected_n_bytes).unwrap();
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    BufWriter::new(universe_history_file)
}

pub fn write_historic_snapshot<T: Write>(
    universe_history_writer: &mut BufWriter<T>,
    universe: &Universe,
    current_time: f64,
    time_step: f64,
) {
    // It can be excessively inefficient to work directly with something that implements Write. For
    // example, every call to write on File results in a system call. A BufWriter keeps an
    // in-memory buffer of data and writes it to an underlying writer in large, infrequent batches.
    //
    // The buffer will be written out when the writer is dropped.

    // 2.- Write accumulative output data to conserve the history of the simulation
    for (current_particle_index, particle) in universe.particles[..universe.n_particles]
        .iter()
        .enumerate()
    {
        // Serialize in chunks of maximum 12 elements or it fails
        let output = (
            current_time, // days
            time_step,    // days
            (particle.id as i32),
            particle.inertial_position.x, // AU
            particle.inertial_position.y,
            particle.inertial_position.z,
            particle.spin.x,
            particle.spin.y,
            particle.spin.z,
            particle.inertial_velocity.x, // AU/days
            particle.inertial_velocity.y,
            particle.inertial_velocity.z,
            particle.mass,   // Msun
            particle.radius, // Rsun
            particle.radius_of_gyration_2,
        );
        let love_number = match &particle.tides.effect {
            TidesEffect::CentralBody(tidal_model) | TidesEffect::OrbitingBody(tidal_model) => {
                match tidal_model {
                    TidalModel::ConstantTimeLag(params) => params.love_number,
                    _ => 0.,
                }
            }
            TidesEffect::Disabled => 0.,
        };
        bincode::serialize_into(&mut *universe_history_writer, &output).unwrap();
        let output = (
            love_number,
            particle.tides.parameters.internal.scaled_dissipation_factor,
            particle.tides.parameters.internal.lag_angle,
            particle.tides.parameters.internal.denergy_dt, // Msun.AU^2.day^-3
            particle.disk.parameters.internal.migration_timescale,
        );
        bincode::serialize_into(&mut *universe_history_writer, &output).unwrap();

        if MIN_ORBITAL_PERIOD_TIME_STEP_RATIO > 0. {
            let reference_particle_index;
            let (
                _semimajor_axis,
                _perihelion_distance,
                _eccentricity,
                _inclination,
                _longitude_of_perihelion,
                _longitude_of_ascending_node,
                _mean_anomaly,
                orbital_period,
            ) = match particle.reference {
                Reference::MostMassiveParticle => {
                    reference_particle_index = universe.hosts.index.most_massive;
                    let reference_particle = &universe.particles[reference_particle_index];
                    calculate_keplerian_orbital_elements(
                        reference_particle.mass_g + particle.mass_g,
                        particle.heliocentric_position,
                        particle.heliocentric_velocity,
                    )
                }
                Reference::Particle(index) => {
                    reference_particle_index = index;
                    let reference_particle = &universe.particles[reference_particle_index];
                    let position = Axes::from(
                        particle.inertial_position.x - reference_particle.inertial_position.x,
                        particle.inertial_position.y - reference_particle.inertial_position.y,
                        particle.inertial_position.z - reference_particle.inertial_position.z,
                    );
                    let velocity = Axes::from(
                        particle.inertial_velocity.x - reference_particle.inertial_velocity.x,
                        particle.inertial_velocity.y - reference_particle.inertial_velocity.y,
                        particle.inertial_velocity.z - reference_particle.inertial_velocity.z,
                    );
                    calculate_keplerian_orbital_elements(
                        reference_particle.mass_g + particle.mass_g,
                        position,
                        velocity,
                    )
                }
            };

            // Control once in a while (when historic point is written) that the
            // time step is small enough to correctly integrate an orbit
            if current_particle_index != reference_particle_index
                && orbital_period <= time_step * MIN_ORBITAL_PERIOD_TIME_STEP_RATIO
            {
                println!("\n");
                panic!(
                    "[PANIC {} UTC] Time step is too large! Particle {} has an orbital period around particle {} of {:0.3} days which is less than the recommended limit ({:0.3} days) based on the current time step ({:0.3} days).",
                    OffsetDateTime::now_utc()
                        .format(
                            &format_description::parse(
                                "[year].[month].[day] [hour]:[minute]:[second]"
                            )
                            .unwrap()
                        )
                        .unwrap(),
                    current_particle_index,
                    reference_particle_index,
                    orbital_period,
                    time_step * MIN_ORBITAL_PERIOD_TIME_STEP_RATIO,
                    time_step
                );
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
//- Restore functions
////////////////////////////////////////////////////////////////////////////////

pub fn restore_snapshot(
    universe_integrator_snapshot_path: &Path,
) -> Result<Box<dyn Integrator>, String> {
    let mut universe_integrator: Box<dyn Integrator>;
    if universe_integrator_snapshot_path.exists() {
        // Open the path in read-only mode only to verify it exists, returns `io::Result<File>`
        if let Err(why) = File::open(universe_integrator_snapshot_path) {
            return Err(format!(
                "Couldn't open {}: {}",
                universe_integrator_snapshot_path.display(),
                why
            ));
        }

        if universe_integrator_snapshot_path.extension().unwrap() == "json" {
            universe_integrator =
                deserialize_json_snapshot(universe_integrator_snapshot_path).unwrap();
        } else {
            universe_integrator =
                deserialize_bin_snapshot(universe_integrator_snapshot_path).unwrap();
        }
        if universe_integrator.get_current_time() == 0. {
            println!(
                "[INFO {} UTC] Created new simulation based on '{}'.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                universe_integrator_snapshot_path.display()
            );
            universe_integrator.initialize_physical_values();
        } else {
            println!(
                "[INFO {} UTC] Restored previous simulation from '{}'.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                universe_integrator_snapshot_path.display()
            );
            let current_time_years = universe_integrator.get_current_time() / 365.25;
            println!(
                "[INFO {} UTC] Continuing from year {:0.0} ({:0.1e}).",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap(),
                current_time_years,
                current_time_years
            );
        }
        Ok(universe_integrator)
    } else {
        Err("File does not exist".to_string())
    }
}

fn deserialize_json_snapshot(snapshot_path: &Path) -> Result<Box<dyn Integrator>, String> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let mut snapshot_file = File::open(snapshot_path).unwrap();
    //// Deserialize using `json::decode`
    let mut json_encoded = String::new();
    if let Err(why) = snapshot_file.read_to_string(&mut json_encoded) {
        return Err(format!("Couldn't read json snapshot file: {why}"));
    }

    let wrapped_universe_integrator: Result<WHFast, serde_json::Error> =
        serde_json::from_str(&json_encoded);
    if let Ok(universe_integrator) = wrapped_universe_integrator {
        println!(
            "[INFO {} UTC] WHFAST Integrator.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        Ok(Box::new(universe_integrator))
    } else {
        let wrapped_universe_integrator: Result<Ias15, serde_json::Error> =
            serde_json::from_str(&json_encoded);
        if let Ok(universe_integrator) = wrapped_universe_integrator {
            println!(
                "[INFO {} UTC] IAS15 Integrator.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            Ok(Box::new(universe_integrator))
        } else {
            let wrapped_universe_integrator: Result<LeapFrog, serde_json::Error> =
                serde_json::from_str(&json_encoded);
            match wrapped_universe_integrator {
                Ok(universe_integrator) => {
                    println!(
                        "[INFO {} UTC] LeapFrog Integrator.",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap()
                    );
                    Ok(Box::new(universe_integrator))
                }
                Err(_) => Err("Unknown integrator!".to_string()),
            }
        }
    }
}

fn deserialize_bin_snapshot(snapshot_path: &Path) -> Result<Box<dyn Integrator>, String> {
    // Open the path in read-only mode, returns `io::Result<File>`
    let snapshot_file = File::open(snapshot_path).unwrap();
    let mut reader = BufReader::new(&snapshot_file);
    let wrapped_universe_integrator: Result<WHFast, bincode::Error> =
        bincode::deserialize_from(&mut reader);
    if let Ok(universe_integrator) = wrapped_universe_integrator {
        println!(
            "[INFO {} UTC] WHFAST Integrator.",
            OffsetDateTime::now_utc()
                .format(
                    &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                        .unwrap()
                )
                .unwrap()
        );
        Ok(Box::new(universe_integrator))
    } else {
        // Re-open file because the previous File/BufReader was already consumed
        let snapshot_file = File::open(snapshot_path).unwrap();
        let mut reader = BufReader::new(&snapshot_file);
        let wrapped_universe_integrator: Result<Ias15, bincode::Error> =
            bincode::deserialize_from(&mut reader);
        if let Ok(universe_integrator) = wrapped_universe_integrator {
            println!(
                "[INFO {} UTC] IAS15 Integrator.",
                OffsetDateTime::now_utc()
                    .format(
                        &format_description::parse("[year].[month].[day] [hour]:[minute]:[second]")
                            .unwrap()
                    )
                    .unwrap()
            );
            Ok(Box::new(universe_integrator))
        } else {
            // Re-open file because the previous File/BufReader was already consumed
            let snapshot_file = File::open(snapshot_path).unwrap();
            let mut reader = BufReader::new(&snapshot_file);
            let wrapped_universe_integrator: Result<LeapFrog, bincode::Error> =
                bincode::deserialize_from(&mut reader);
            match wrapped_universe_integrator {
                Ok(universe_integrator) => {
                    println!(
                        "[INFO {} UTC] LeapFrog Integrator.",
                        OffsetDateTime::now_utc()
                            .format(
                                &format_description::parse(
                                    "[year].[month].[day] [hour]:[minute]:[second]"
                                )
                                .unwrap()
                            )
                            .unwrap()
                    );
                    Ok(Box::new(universe_integrator))
                }
                Err(_) => Err("Unknown integrator!".to_string()),
            }
        }
    }
}
