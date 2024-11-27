pub mod stars;
pub mod planets;
pub mod universe;

use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufRead};

#[allow(dead_code)]
pub fn simulation_properties() -> (f64, f64, f64, f64, f64) {
    let time_step: f64 = 0.08; // in days
    let time_limit: f64 = time_step*200.; // days
    let initial_time: f64 = 1.6e8*365.25; // time [days] where simulation starts
    let historic_snapshot_period: f64 = 100.*365.25; // days
    let recovery_snapshot_period: f64 = 10.*historic_snapshot_period; // days
    (time_step, time_limit, initial_time, historic_snapshot_period, recovery_snapshot_period)
}

#[allow(dead_code)]
pub fn get_data_dirname(test_name: &String) -> (String, String) {
    let rust_data_dirname = format!("tests/data/{0}/", test_name);
    let python_data_dirname = format!("posidonius/tests/data/{0}/", test_name);
    (rust_data_dirname, python_data_dirname)
}

#[allow(dead_code)]
pub fn load_kaula_parameters(file_path: &str) -> Result<posidonius::KaulaParameters, Box<dyn Error>> {
    // Open the file
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // Temporary vectors to hold the columns before assigning them to the struct
    let mut w_lm: Vec<f64> = Vec::new();
    let mut im_k2: Vec<f64> = Vec::new();
    let mut re_k2: Vec<f64> = Vec::new();

    for line in reader.lines() {
        let line = line?;

        // Skip lines that are comments (starting with '#')
        if line.starts_with('#') {
            continue;
        }

        // Use split_whitespace() to handle any number of spaces between columns
        let parts: Vec<&str> = line.split_whitespace().collect();

        // Ensure we have at least 3 columns
        if parts.len() >= 3 {
            w_lm.push(parts[0].parse()?);    // Column 0
            im_k2.push(parts[1].parse()?);  // Column 1
            re_k2.push(parts[2].parse()?);  // Column 2
        }
    }

    // NewFeatures: Adding zeros to the array is not exactly a good practice, weird result was encountered
    // when comparing the scientific test cases, in which I extened the array with the last values.
    // Perhaps the result coming from the unsorted arrays followed by the binary search in kaula.rs.
    // This need to be modified if one we decided to go for a dynamical sizes for love number spectra.
    let expected_size = 32 * 32;

    // Get the number of available data points
    let num_datapoints = w_lm.len();

    if num_datapoints > expected_size {
        panic!(
            "Data size exceeds expected size: num_datapoints ({}) > expected_size ({})",
            num_datapoints, expected_size
        );
    }

    // Extend w_lm, im_k2, and re_k2 to expected_size by repeating the last element if needed
    if num_datapoints < expected_size {
        let last_w_lm = *w_lm.last().unwrap_or(&0.0);
        let last_im_k2 = *im_k2.last().unwrap_or(&0.0);
        let last_re_k2 = *re_k2.last().unwrap_or(&0.0);

        w_lm.resize(expected_size, last_w_lm);
        im_k2.resize(expected_size, last_im_k2);
        re_k2.resize(expected_size, last_re_k2);
    }

    // Create a zeroed out array for the Kaula parameters (size: 32*32 = 1024)
    let mut love_number_excitation_frequency = [0.0; 32 * 32];
    let mut real_part_love_number = [0.0; 32 * 32];
    let mut imaginary_part_love_number = [0.0; 32 * 32];

    // Copy data from the vectors into the fixed-size arrays
    love_number_excitation_frequency.copy_from_slice(&w_lm[..expected_size]);
    real_part_love_number.copy_from_slice(&re_k2[..expected_size]);
    imaginary_part_love_number.copy_from_slice(&im_k2[..expected_size]);

    // Build and return the KaulaParameters struct
    Ok(posidonius::KaulaParameters {
        love_number_excitation_frequency,
        real_part_love_number,
        imaginary_part_love_number,
        num_datapoints: num_datapoints as f64,
        // below two are placeholder for the moment
        stellar_tide: -1, // Specified the type of tides
        spectrum_spin_rate: 0.0, // Specified the initial spin rate of the star (stellar tide specific)
        polynomials: posidonius::Polynomials::new(),
        kaula_tidal_force: posidonius::Axes { x: 0.0, y: 0.0, z: 0.0 },
    })
}
