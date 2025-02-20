use super::super::constants::SUN_DYN_FREQ_2;
use super::super::constants::{M2AU, R_SUN};
use super::super::tools::linear_interpolation;
use super::super::{Particle, TidalModel, TidesEffect};
use csv;
use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum EvolutionType {
    GalletBolmont2017(f64), // SolarLike EvolvingDissipation Evolving dissipation
    BolmontMathis2016(f64), // SolarLike EvolvingDissipation Evolving dissipation
    Baraffe2015(f64),       // NEW
    Leconte2011(f64),       // BrownDwarf
    Baraffe1998(f64),       // M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
    LeconteChabrier2013(bool), // Jupiter with/without dissipation of dynamical tides
    NonEvolving,
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct Evolver {
    pub evolution: EvolutionType,
    pub time: Vec<f64>,
    pub radius: Vec<f64>,
    pub radius_of_gyration_2: Vec<f64>,
    pub love_number: Vec<f64>,
    pub inverse_tidal_q_factor: Vec<f64>, // Bolmont & Mathis 2016
    left_index: usize,
}

// NOTE: This is a big optimization to reduce the cost of cloning
//       Cloned particles are downgraded to a NonEvolving status
//       to avoid cloning model's data which can be big and expensive.
impl Clone for Evolver {
    fn clone(&self) -> Self {
        Evolver {
            evolution: EvolutionType::NonEvolving,
            time: vec![],
            radius: vec![],
            radius_of_gyration_2: vec![],
            love_number: vec![],
            inverse_tidal_q_factor: vec![],
            left_index: 0,
        }
    }
}

impl Evolver {
    pub fn new(evolution: EvolutionType, initial_time: f64, time_limit: f64) -> Evolver {
        let mut time = vec![];
        let mut radius = vec![];
        let mut radius_of_gyration_2 = vec![];
        let mut love_number = vec![];
        let mut inverse_tidal_q_factor = vec![];

        let filename = match evolution {
            EvolutionType::GalletBolmont2017(mass) => {
                if (0.299..=0.301).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_03_Z_0134.dat")
                } else if (0.399..=0.401).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_04_Z_0134.dat")
                } else if (0.599..=0.601).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_06_Z_0134.dat")
                } else if (0.699..=0.701).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_07_Z_0134.dat")
                } else if (0.799..=0.801).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_08_Z_0134.dat")
                } else if (0.899..=0.901).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_09_Z_0134.dat")
                } else if (0.999..=1.001).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_10_Z_0134.dat")
                } else if (1.099..=1.101).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_11_Z_0134.dat")
                } else if (1.199..=1.201).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_12_Z_0134.dat")
                } else if (1.399..=1.401).contains(&mass) {
                    String::from("input/Gallet_Bolmont_2017/M_14_Z_0134.dat")
                } else {
                    panic!(
                        "The evolution type Gallet_Bolmont_2017 does not support a mass of {mass} Msun!"
                    );
                }
            }
            EvolutionType::Baraffe2015(mass) => {
                if (0.0099..=0.0101).contains(&mass) {
                    String::from("input/Baraffe_2015/0010_Msun.dat")
                } else if (0.0149..=0.0151).contains(&mass) {
                    String::from("input/Baraffe_2015/0015_Msun.dat")
                } else if (0.0199..=0.0201).contains(&mass) {
                    String::from("input/Baraffe_2015/0020_Msun.dat")
                } else if (0.0299..=0.0301).contains(&mass) {
                    String::from("input/Baraffe_2015/0030_Msun.dat")
                } else if (0.0399..=0.0401).contains(&mass) {
                    String::from("input/Baraffe_2015/0040_Msun.dat")
                } else if (0.0499..=0.0501).contains(&mass) {
                    String::from("input/Baraffe_2015/0050_Msun.dat")
                } else if (0.0599..=0.0601).contains(&mass) {
                    String::from("input/Baraffe_2015/0060_Msun.dat")
                } else if (0.0699..=0.0701).contains(&mass) {
                    String::from("input/Baraffe_2015/0070_Msun.dat")
                } else if (0.0719..=0.0721).contains(&mass) {
                    String::from("input/Baraffe_2015/0072_Msun.dat")
                } else if (0.0749..=0.0751).contains(&mass) {
                    String::from("input/Baraffe_2015/0075_Msun.dat")
                } else if (0.0799..=0.0801).contains(&mass) {
                    String::from("input/Baraffe_2015/0080_Msun.dat")
                } else if abs!(mass - 0.09) < 1e-7 {
                    String::from("input/Baraffe_2015/0090_Msun.dat")
                } else if abs!(mass - 0.11) < 1e-7 {
                    String::from("input/Baraffe_2015/0110_Msun.dat")
                } else if abs!(mass - 0.13) < 1e-7 {
                    String::from("input/Baraffe_2015/0130_Msun.dat")
                } else if abs!(mass - 0.15) < 1e-7 {
                    String::from("input/Baraffe_2015/0150_Msun.dat")
                } else if abs!(mass - 0.17) < 1e-7 {
                    String::from("input/Baraffe_2015/0170_Msun.dat")
                } else if abs!(mass - 0.20) < 1e-7 {
                    String::from("input/Baraffe_2015/0200_Msun.dat")
                } else if abs!(mass - 0.30) < 1e-7 {
                    String::from("input/Baraffe_2015/0300_Msun.dat")
                } else if abs!(mass - 0.40) < 1e-7 {
                    String::from("input/Baraffe_2015/0400_Msun.dat")
                } else if abs!(mass - 0.50) < 1e-7 {
                    String::from("input/Baraffe_2015/0500_Msun.dat")
                } else if abs!(mass - 0.60) < 1e-7 {
                    String::from("input/Baraffe_2015/0600_Msun.dat")
                } else if abs!(mass - 0.70) < 1e-7 {
                    String::from("input/Baraffe_2015/0700_Msun.dat")
                } else if abs!(mass - 0.80) < 1e-7 {
                    String::from("input/Baraffe_2015/0800_Msun.dat")
                } else if abs!(mass - 0.90) < 1e-7 {
                    String::from("input/Baraffe_2015/0900_Msun.dat")
                } else if abs!(mass - 1.00) < 1e-7 {
                    String::from("input/Baraffe_2015/1000_Msun.dat")
                } else if abs!(mass - 1.10) < 1e-7 {
                    String::from("input/Baraffe_2015/1100_Msun.dat")
                } else if abs!(mass - 1.20) < 1e-7 {
                    String::from("input/Baraffe_2015/1200_Msun.dat")
                } else if abs!(mass - 1.30) < 1e-7 {
                    String::from("input/Baraffe_2015/1300_Msun.dat")
                } else if abs!(mass - 1.40) < 1e-7 {
                    String::from("input/Baraffe_2015/1400_Msun.dat")
                } else {
                    panic!(
                        "The evolution type Baraffe2015 does not support a mass of {mass} Msun!",
                    );
                }
            }
            EvolutionType::Leconte2011(mass) => {
                if (0.0099..=0.0101).contains(&mass) {
                    String::from("input/Leconte_2011/mass_10.0000.dat")
                } else if (0.0119..=0.0121).contains(&mass) {
                    String::from("input/Leconte_2011/mass_12.0000.dat")
                } else if (0.0149..=0.0151).contains(&mass) {
                    String::from("input/Leconte_2011/mass_15.0000.dat")
                } else if (0.0199..=0.0201).contains(&mass) {
                    String::from("input/Leconte_2011/mass_20.0000.dat")
                } else if (0.0299..=0.0301).contains(&mass) {
                    String::from("input/Leconte_2011/mass_30.0000.dat")
                } else if (0.0399..=0.0401).contains(&mass) {
                    String::from("input/Leconte_2011/mass_40.0000.dat")
                } else if (0.0499..=0.0501).contains(&mass) {
                    String::from("input/Leconte_2011/mass_50.0000.dat")
                } else if (0.0599..=0.0601).contains(&mass) {
                    String::from("input/Leconte_2011/mass_60.0000.dat")
                } else if (0.0699..=0.0701).contains(&mass) {
                    String::from("input/Leconte_2011/mass_70.0000.dat")
                } else if (0.0719..=0.0721).contains(&mass) {
                    String::from("input/Leconte_2011/mass_72.0000.dat")
                } else if (0.0749..=0.0751).contains(&mass) {
                    String::from("input/Leconte_2011/mass_75.0000.dat")
                } else if (0.0799..=0.0801).contains(&mass) {
                    String::from("input/Leconte_2011/mass_80.0000.dat")
                } else {
                    panic!(
                        "The evolution type Leconte2011 does not support a mass of {mass} Msun!",
                    );
                }
            }
            EvolutionType::Baraffe1998(mass) => {
                if abs!(mass - 0.10) <= 1.0e-7 {
                    String::from("input/Baraffe_1998/01Msun.dat")
                } else if abs!(mass - 1.0) <= 1.0e-7 {
                    String::from("input/Baraffe_1998/SRad_Spli_M-1_0000.dat")
                } else {
                    panic!(
                        "The evolution type Baraffe1998 does not support a mass of {mass} Msun!"
                    );
                }
            }

            EvolutionType::BolmontMathis2016(mass) => {
                if (0.399..=0.401).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L04Z02r.dat")
                } else if (0.499..=0.501).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L05Z02r.dat")
                } else if (0.599..=0.601).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L06Z02r.dat")
                } else if (0.699..=0.701).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L07Z02r.dat")
                } else if (0.799..=0.801).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L08Z02r.dat")
                } else if (0.899..=0.901).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L09Z02r.dat")
                } else if (0.999..=1.001).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L10Z02r.dat")
                } else if (1.099..=1.101).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L11Z02r.dat")
                } else if (1.199..=1.201).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L12Z02r.dat")
                } else if (1.299..=1.301).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L13Z02r.dat")
                } else if (1.399..=1.401).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L14Z02r.dat")
                } else if (1.499..=1.501).contains(&mass) {
                    String::from("input/Bolmont_Mathis_2016/L15Z02r.dat")
                } else {
                    panic!(
                        "The evolution type BolmontMathis2016 does not support a mass of {mass} Msun!",
                    );
                }
            }
            EvolutionType::LeconteChabrier2013(_) => {
                String::from("input/Leconte_Chabrier_2013/Jupiter.dat")
            }
            EvolutionType::NonEvolving => String::from("input/empty.dat"),
        };
        //println!("Filename {}", filename);

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b' ')
            .flexible(true)
            .from_path(&filename)
            .unwrap();
        for (i, row) in rdr.records().map(|r| r.unwrap()).enumerate() {
            let raw_time = row[0].parse::<f64>().unwrap_or(-1.);
            if i == 0 && raw_time < 0. {
                // First line, probably column headers
                continue;
            }
            let raw_radius = match evolution {
                EvolutionType::Baraffe2015(_) => row[2].parse::<f64>().unwrap(),
                EvolutionType::GalletBolmont2017(_) => row[3].parse::<f64>().unwrap(),
                _ => row[1].parse::<f64>().unwrap(),
            };
            // All types have time and radius
            let (current_time, current_radius) = match evolution {
                EvolutionType::Leconte2011(_) => (raw_time, raw_radius * R_SUN),
                EvolutionType::Baraffe1998(mass) => {
                    if abs!(mass - 0.10) <= 1.0e-7 {
                        (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                    } else if abs!(mass - 1.0) <= 1.0e-7 {
                        (raw_time * 365.25 - initial_time, raw_radius * M2AU)
                    } else {
                        (0., 0.)
                    }
                }
                EvolutionType::Baraffe2015(_) | EvolutionType::BolmontMathis2016(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * R_SUN)
                }
                EvolutionType::GalletBolmont2017(_) => (
                    10_f64.powf(raw_time) * 365.25 - initial_time,
                    raw_radius * R_SUN,
                ),
                EvolutionType::LeconteChabrier2013(_) => {
                    (raw_time * 365.25 - initial_time, raw_radius * M2AU)
                }
                EvolutionType::NonEvolving => (0., 0.),
            };
            // Fields that only some types have...
            let current_radius_of_gyration_2 = match evolution {
                EvolutionType::LeconteChabrier2013(_) | EvolutionType::Baraffe2015(_) => {
                    row[3].parse::<f64>().unwrap()
                }
                _ => 0.,
            };
            let current_love_number = match evolution {
                EvolutionType::LeconteChabrier2013(_) => row[2].parse::<f64>().unwrap(),
                _ => 0.,
            };
            let current_inverse_tidal_q_factor = match evolution {
                EvolutionType::BolmontMathis2016(_) => row[2].parse::<f64>().unwrap(),
                EvolutionType::GalletBolmont2017(_) => {
                    let raw_inverse_tidal_q_factor = row[10].parse::<f64>().unwrap();
                    1. / 10_f64.powf(raw_inverse_tidal_q_factor)
                }
                EvolutionType::LeconteChabrier2013(true) => row[5].parse::<f64>().unwrap(),
                _ => 0.,
            };

            // Only save values needed for the evolution model to save memory
            // and computation time when cloning
            match evolution {
                EvolutionType::GalletBolmont2017(_) | EvolutionType::BolmontMathis2016(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                }
                EvolutionType::Baraffe2015(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                }
                EvolutionType::Leconte2011(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    // Read radius of gyration from an auxiliary file later on
                    //radius_of_gyration_2.push(current_radius_of_gyration_2);
                }
                EvolutionType::Baraffe1998(_) => {
                    time.push(current_time);
                    radius.push(current_radius);
                }
                EvolutionType::LeconteChabrier2013(dissipation_of_dynamical_tides) => {
                    time.push(current_time);
                    radius.push(current_radius);
                    radius_of_gyration_2.push(current_radius_of_gyration_2);
                    love_number.push(current_love_number);
                    if dissipation_of_dynamical_tides {
                        inverse_tidal_q_factor.push(current_inverse_tidal_q_factor);
                    }
                }
                EvolutionType::NonEvolving => {}
            }
        }

        // Leconte2011 have a separate file for radius of gyration with a different time sampling:
        if let EvolutionType::Leconte2011(mass) = evolution {
            // Read separate file
            let mut aux_time = vec![];
            let mut aux_radius_of_gyration_2 = vec![];
            // The column to read from the auxiliary file depends on the mass
            let aux_column = {
                if (0.0099..=0.0101).contains(&mass) {
                    1
                } else if (0.0119..=0.0121).contains(&mass) {
                    2
                } else if (0.0149..=0.0151).contains(&mass) {
                    3
                } else if (0.0199..=0.0201).contains(&mass) {
                    4
                } else if (0.0299..=0.0301).contains(&mass) {
                    5
                } else if (0.0399..=0.0401).contains(&mass) {
                    6
                } else if (0.0499..=0.0501).contains(&mass) {
                    7
                } else if (0.0599..=0.0601).contains(&mass) {
                    8
                } else if (0.0699..=0.0701).contains(&mass) {
                    9
                } else if (0.0719..=0.0721).contains(&mass) {
                    10
                } else if (0.0749..=0.0751).contains(&mass) {
                    11
                } else if (0.0799..=0.0801).contains(&mass) {
                    12
                } else {
                    panic!(
                        "The evolution type Leconte2011 does not support a mass of {mass} Msun!",
                    );
                }
            };
            let mut rdr = csv::ReaderBuilder::new()
                .has_headers(false)
                .delimiter(b' ')
                .flexible(true)
                .from_path("input/Leconte_2011/rg2BD.dat")
                .unwrap();
            for row in rdr.records().map(|r| r.unwrap()) {
                let raw_time = row[0].parse::<f64>().unwrap();
                let raw_radius_of_gyration_2 = row[aux_column].parse::<f64>().unwrap();
                aux_time.push(raw_time);
                aux_radius_of_gyration_2.push(raw_radius_of_gyration_2);
            }
            // The time range from the main file and the auxiliary do not match
            // the main one is bigger, so we restrain to the smallest one (the aux one)
            if let Some(lower_limit) = time.iter().position(|&x| x >= *aux_time.first().unwrap()) {
                time.drain(0..lower_limit);
                radius.drain(0..lower_limit);
            }
            if let Some(upper_limit) = time.iter().position(|&x| x > *aux_time.last().unwrap()) {
                time.drain(upper_limit..);
                radius.drain(upper_limit..);
            }

            // Interpolate the radius of gyration to match the same time sampling shared by other parameters (e.g., radius)
            let precision = 7;
            for current_time in &time {
                let (current_radius_of_gyration_2, _) =
                    linear_interpolation(*current_time, &aux_time, &aux_radius_of_gyration_2);
                radius_of_gyration_2.push(math::round::half_up(
                    current_radius_of_gyration_2,
                    precision,
                ));
                //println!("Rg2 {} {}", current_time, current_radius_of_gyration_2);
            }

            // Convert time to be homogenous with the rest of models
            time = time.iter().map(|&t| t * 365.25 - initial_time).collect();
        }

        assert!(
            time.is_empty() || time[0] <= 0.,
            "Your initial time ({} days) is smaller than the minimum allowed age of the star ({} days)",
            initial_time,
            time[0] + initial_time
        );
        assert!(
            time.is_empty() || time[time.len() - 1] >= time_limit,
            "Your time limit ({} days) is greater than the maximum allowed age of the star ({} days)",
            time_limit,
            time[time.len() - 1]
        );

        Evolver {
            evolution,
            time,
            radius,
            radius_of_gyration_2,
            love_number,
            inverse_tidal_q_factor,
            left_index: 0,
        }
    }

    // OPTIMIZATION: Skip first N elements which belong to the past
    // but still go through the last n elements just in case the integrator
    // is travelling to the past (e.g., current time step did not converge)
    fn idx(&self) -> usize {
        let security_shift = 10;
        if self.left_index > security_shift {
            self.left_index - security_shift
        } else {
            0
        }
    }

    pub fn radius(&mut self, current_time: f64, current_radius: f64) -> f64 {
        let (new_radius, left_index) = match self.evolution {
            EvolutionType::NonEvolving => (current_radius, 0),
            _ => linear_interpolation(
                current_time,
                &self.time[self.idx()..],
                &self.radius[self.idx()..],
            ),
        };
        self.left_index += left_index;
        new_radius
    }

    pub fn radius_of_gyration_2(
        &mut self,
        current_time: f64,
        current_radius_of_gyration_2: f64,
    ) -> f64 {
        match self.evolution {
            EvolutionType::Baraffe2015(_)
            | EvolutionType::Leconte2011(_)
            | EvolutionType::LeconteChabrier2013(_) => {
                let (new_radius_of_gyration_2, left_index) = linear_interpolation(
                    current_time,
                    &self.time[self.idx()..],
                    &self.radius_of_gyration_2[self.idx()..],
                );
                self.left_index += left_index;

                new_radius_of_gyration_2
            }
            _ => current_radius_of_gyration_2,
        }
    }

    pub fn love_number(&mut self, current_time: f64, current_love_number: f64) -> f64 {
        let (new_love_number, left_index) = match self.evolution {
            EvolutionType::LeconteChabrier2013(_) => linear_interpolation(
                current_time,
                &self.time[self.idx()..],
                &self.love_number[self.idx()..],
            ),
            _ => (current_love_number, 0),
        };
        self.left_index += left_index;
        new_love_number
    }

    pub fn inverse_tidal_q_factor(
        &mut self,
        current_time: f64,
        current_inverse_tidal_q_factor: f64,
    ) -> f64 {
        // The planet excite/induce waves in the convective region of the star
        // and the Q factor is a quality factor that describes that loss of energy.
        //
        // The inverse_tidal_q_factor was derived from analytical formulas using
        // parameters from stellar models. This factor is used to derive the
        // dissipation factor of the star, which depends on the excitation
        // frequency of the planet (the frequency at what we would see the planet
        // if we were sitting at the surface of the star).
        //
        //  - If the star does not rotate, the excitation frequency is the orbital period (i.e.,
        //  frequency)
        //  - If the planet is very far away, the excitation frequency is the spin
        //  - Planets are treated as source points and their rotation is not important here
        //
        // The tidal theory (and this code) assumes circular orbits because they are easier.
        // Excentric orbits needs more than one frequency to be described and it will be
        // included in future versions of this code.
        let (new_inverse_tidal_q_factor, left_index) = match self.evolution {
            EvolutionType::BolmontMathis2016(_)
            | EvolutionType::GalletBolmont2017(_)
            | EvolutionType::LeconteChabrier2013(true) => linear_interpolation(
                current_time,
                &self.time[self.idx()..],
                &self.inverse_tidal_q_factor[self.idx()..],
            ),
            _ => (current_inverse_tidal_q_factor, 0),
        };
        self.left_index += left_index;
        new_inverse_tidal_q_factor
    }
}

pub fn calculate_particles_non_spin_dependent_evolving_quantities(
    current_time: f64,
    particles: &mut [Particle],
    particles_evolvers: &mut [Evolver],
) {
    for (particle, evolver) in particles.iter_mut().zip(particles_evolvers.iter_mut()) {
        ////////////////////////////////////////////////////////////////////
        // Radius and radius of gyration 2
        ////////////////////////////////////////////////////////////////////
        // If particle radius/radius_of_gyration_2 evolves
        // - Compute the moment of inertia
        let new_radius = evolver.radius(current_time, particle.radius);
        let new_radius_of_gyration_2 =
            evolver.radius_of_gyration_2(current_time, particle.radius_of_gyration_2);
        if new_radius != particle.radius
            || new_radius_of_gyration_2 != particle.radius_of_gyration_2
        {
            particle.radius = new_radius;
            particle.radius_of_gyration_2 = new_radius_of_gyration_2;
            particle.moment_of_inertia =
                particle.mass * particle.radius_of_gyration_2 * particle.radius.powi(2);
        }

        ////////////////////////////////////////////////////////////////////
        // Love number
        ////////////////////////////////////////////////////////////////////
        match &mut particle.tides.effect {
            &mut (TidesEffect::CentralBody(mut tidal_model)
            | TidesEffect::OrbitingBody(mut tidal_model)) => {
                if let &mut TidalModel::ConstantTimeLag(mut params) = &mut tidal_model {
                    params.love_number = evolver.love_number(current_time, params.love_number);
                }
            }
            TidesEffect::Disabled => {}
        }

        //println!("[{}] Evolve Radius {:e} Gyration {:e} Love {:e}", current_time, particle.radius, particle.radius_of_gyration_2, particle.love_number);
    }
}

pub fn calculate_particles_spin_dependent_evolving_quantities(
    current_time: f64,
    particles: &mut [Particle],
    particles_evolvers: &mut [Evolver],
) {
    for (particle, evolver) in particles.iter_mut().zip(particles_evolvers.iter_mut()) {
        ////////////////////////////////////////////////////////////////////
        // Lag angle
        ////////////////////////////////////////////////////////////////////
        particle.tides.parameters.internal.lag_angle = match evolver.evolution {
            EvolutionType::BolmontMathis2016(_)
            | EvolutionType::GalletBolmont2017(_)
            | EvolutionType::LeconteChabrier2013(true) => {
                let inverse_tidal_q_factor = evolver.inverse_tidal_q_factor(current_time, 0.);
                let epsilon_squared = particle.norm_spin_vector_2 / SUN_DYN_FREQ_2;
                // Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
                // but as for sigma it is necessary to divide by k2s, we do not divide here

                3.0 * epsilon_squared * inverse_tidal_q_factor / 4.0
            }
            _ => 0.,
        };
        //println!("[{}] Evolve Lag {:e}", current_time, particle.lag_angle);
    }
}
