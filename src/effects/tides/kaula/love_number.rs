use serde::{Deserialize, Serialize};
use serde_big_array::BigArray;

use super::select_eccentricty_order_q;
use crate::constants::MAX_PARTICLES;

// Relative change of the spin rate or of the orbital frequency below which the cached love
// numbers are not recomputed (the tidal frequencies moved by less than this fraction; the
// spectrum is interpolated on a grid ~1e-6 rad/s wide, so the k2 values do not change).
const LOVE_NUMBER_REFRESH_TOLERANCE: f64 = 1e-8;

// Tidal frequency, real and imaginary love numbers calculated each timestep.
// One cache per perturbing body (the star's spectrum is shared by all its planets, which have
// different orbital frequencies); `valid` and the `last_*` values gate the refresh.
#[derive(Copy, PartialEq, Debug, Clone)]
struct Cache {
    freq: [f64; 145],
    real: [f64; 145],
    imag: [f64; 145],
    valid: bool,
    last_spin_rate: f64,
    last_orbital_frequency: f64,
    last_q_range: (usize, usize),
    last_full: bool,
}

impl Cache {
    const fn empty() -> Self {
        Cache {
            freq: [0.; 145],
            real: [0.; 145],
            imag: [0.; 145],
            valid: false,
            last_spin_rate: 0.,
            last_orbital_frequency: 0.,
            last_q_range: (0, 0),
            last_full: false,
        }
    }

    fn up_to_date(&self, spin_rate: f64, orbital_frequency: f64, q_range: (usize, usize), full: bool) -> bool {
        self.valid
            && self.last_q_range == q_range
            && (self.last_full || !full)
            && (spin_rate - self.last_spin_rate).abs() <= LOVE_NUMBER_REFRESH_TOLERANCE * spin_rate.abs()
            && (orbital_frequency - self.last_orbital_frequency).abs()
                <= LOVE_NUMBER_REFRESH_TOLERANCE * orbital_frequency.abs()
    }
}

impl Cache {
    // Maps a set of 3D coordinates into a 1D index.
    fn map_3d_to_1d(m: usize, m_max: usize, p: usize, p_max: usize, q: usize) -> usize {
        (q * m_max * p_max) + (p * m_max) + m
    }

    /// Fetches the real love number from the cache for index of tuple (m, p, q).
    fn real(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        // m_max and p_max are both 3
        let index = Self::map_3d_to_1d(m, 3, p, 3, q);
        assert!(index < self.real.len());
        self.real[index]
    }

    /// Fetches the imaginary love number from the cache for index of tuple (m, p, q).
    fn imaginary(&self, m: usize, p: usize, q: usize) -> f64 {
        // The cached data is stored in a 1D array, so the 3D coordinates are mapped to the 1D index.
        // m_max and p_max are both 3
        let index = Self::map_3d_to_1d(m, 3, p, 3, q);
        assert!(index < self.imag.len());
        self.imag[index]
    }

    fn set_all(&mut self, m: usize, p: usize, q: usize, frequency: f64, real: f64, imaginary: f64) {
        let index = Self::map_3d_to_1d(m, 3, p, 3, q);
        assert!(index < self.freq.len());
        self.freq[index] = frequency;
        assert!(index < self.real.len());
        self.real[index] = real;
        assert!(index < self.imag.len());
        self.imag[index] = imaginary;
    }
}

#[derive(Serialize, Copy, Deserialize, PartialEq, Debug, Clone)]
#[serde(default, deny_unknown_fields)]
pub struct LoveNumber {
    #[serde(with = "BigArray")]
    spectrum_excitation_frequency: [f64; 1024],
    #[serde(with = "BigArray")]
    spectrum_real_part: [f64; 1024],
    #[serde(with = "BigArray")]
    spectrum_imaginary_part: [f64; 1024],

    // If we are simulating stellar tide, then we need the value of the spectrum spin rate,
    // which was an input value of stellar spin rate from the computation of the spectrum.
    stellar_spectrum_spin_rate: Option<f64>, // Specified the initial spin rate of the star, or `None` if planetary tide.

    // Caches of the tidal frequencies and love numbers, one per perturbing body (indexed by
    // the perturber's particle id); `current` selects the one read by `real`/`imaginary`.
    #[serde(skip)]
    caches: [Cache; MAX_PARTICLES],
    #[serde(skip)]
    current: usize,
}

impl Default for LoveNumber {
    fn default() -> Self {
        Self {
            spectrum_excitation_frequency: [0.; 1024],
            spectrum_real_part: [0.; 1024],
            spectrum_imaginary_part: [0.; 1024],
            stellar_spectrum_spin_rate: None,
            caches: [Cache::empty(); MAX_PARTICLES],
            current: 0,
        }
    }
}

impl LoveNumber {
    pub fn new_from(
        spectrum_excitation_frequency: [f64; 1024],
        spectrum_real_part: [f64; 1024],
        spectrum_imaginary_part: [f64; 1024],
        stellar_spectrum_spin_rate: Option<f64>,
    ) -> Self {
        Self {
            spectrum_excitation_frequency,
            spectrum_real_part,
            spectrum_imaginary_part,
            stellar_spectrum_spin_rate,
            ..Default::default()
        }
    }

    /// Fetches the real love number from the cache for index of tuple (m, p, q).
    pub(crate) fn real(&self, m: usize, p: usize, q: usize) -> f64 {
        self.caches[self.current].real(m, p, q)
    }

    /// Fetches the imaginary love number from the cache for index of tuple (m, p, q).
    pub(crate) fn imaginary(&self, m: usize, p: usize, q: usize) -> f64 {
        self.caches[self.current].imaginary(m, p, q)
    }

    /// `perturber_id`: particle id of the perturbing body, selecting the cache slot.
    pub fn refresh_cache_full(
        &mut self,
        central_body: bool,
        spin_rate: f64,
        orbital_frequency: f64,
        eccentricity: f64,
        perturber_id: usize,
    ) {
        self.current = perturber_id;
        let q_range = select_eccentricty_order_q(eccentricity);
        if self.caches[self.current].up_to_date(spin_rate, orbital_frequency, q_range, true) {
            return;
        }
        // refresh m = 0..3, p = 0..3, q bounded by eccentricty order
        let (m_min, m_max) = (0, 3);
        let (p_min, p_max) = (0, 3);
        self.refresh_cache(
            central_body,
            spin_rate,
            orbital_frequency,
            eccentricity,
            m_min,
            m_max,
            p_min,
            p_max,
        );
        self.caches[self.current] = Cache {
            valid: true,
            last_spin_rate: spin_rate,
            last_orbital_frequency: orbital_frequency,
            last_q_range: q_range,
            last_full: true,
            ..self.caches[self.current]
        };
    }
    /// `perturber_id`: particle id of the perturbing body, selecting the cache slot.
    pub fn refresh_cache_partial(
        &mut self,
        central_body: bool,
        spin_rate: f64,
        orbital_frequency: f64,
        eccentricity: f64,
        perturber_id: usize,
    ) {
        self.current = perturber_id;
        let q_range = select_eccentricty_order_q(eccentricity);
        if self.caches[self.current].up_to_date(spin_rate, orbital_frequency, q_range, false) {
            return;
        }
        // refresh m = 0, p = 1, q = bounded by eccentricty order
        let (m_min, m_max) = (0, 1);
        let (p_min, p_max) = (1, 2);
        self.refresh_cache(
            central_body,
            spin_rate,
            orbital_frequency,
            eccentricity,
            m_min,
            m_max,
            p_min,
            p_max,
        );

        // refresh m = 2, p = 0, q = bounded by eccentricty order
        let (m_min, m_max) = (2, 3);
        let (p_min, p_max) = (0, 1);
        self.refresh_cache(
            central_body,
            spin_rate,
            orbital_frequency,
            eccentricity,
            m_min,
            m_max,
            p_min,
            p_max,
        );
        self.caches[self.current] = Cache {
            valid: true,
            last_spin_rate: spin_rate,
            last_orbital_frequency: orbital_frequency,
            last_q_range: q_range,
            last_full: false,
            ..self.caches[self.current]
        };
    }
    /// Recomputes all the love number values.
    // Called at each time step to cache love numbers for that iteration, to prevent duplicate calculations.
    fn refresh_cache(
        &mut self,
        central_body: bool,
        spin_rate: f64,
        orbital_frequency: f64,
        eccentricity: f64,
        m_min: usize,
        m_max: usize,
        p_min: usize,
        p_max: usize,
    ) {
        let (q_min, q_max) = select_eccentricty_order_q(eccentricity);
        for q in (0..14).skip(q_min).take(q_max) {
            for p in p_min..p_max {
                for m in m_min..m_max {
                    let wk2 = Self::calculate_tidal_excitation_frequency_mode_sigma_2mpq(
                        m,
                        p,
                        q,
                        spin_rate,
                        orbital_frequency,
                    );

                    // TODO Don't update the love numbers cache if the change in wk2 is minimal
                    let (real_k2, imaginary_k2) =
                        self.get_love_numbers(central_body, spin_rate, wk2);

                    self.caches[self.current].set_all(m, p, q, wk2, real_k2, imaginary_k2);
                }
            }
        }
    }

    fn calculate_tidal_excitation_frequency_mode_sigma_2mpq(
        m: usize,
        p: usize,
        q: usize,
        spin: f64,
        orbital_frequency: f64,
    ) -> f64 {
        let (m, p, q) = (m as i32, p as i32, q as i32 - 7);
        f64!(2 - 2 * p + q) * orbital_frequency - f64!(m) * spin
    }

    fn get_love_numbers(&self, central_body: bool, spin_rate: f64, mut wk2: f64) -> (f64, f64) {
        // Only for kaula stellar tide with a fixed spectrum
        // Taking into account of the real time evolution of the stellar spin rate
        // by modifying the tidal frequency, and rescale the love number (below)
        // instead of modifying the spectrum every time step.
        if let Some(spectrum_spin_rate) = self.stellar_spectrum_spin_rate {
            // Days
            wk2 *= spectrum_spin_rate / spin_rate;
        }

        // Planets have symmetric tidal response, stars do not.
        let parity = !central_body & (wk2 < 0.0);
        if parity {
            wk2 = abs!(wk2);
        }

        let (mut real_k2, mut imaginary_k2) = self.interpolate_love_numbers(wk2);

        real_k2 *= -1.;

        // Reverse the parity inversion
        if parity {
            imaginary_k2 *= -1.;
        }

        // Rescale the love number
        if let Some(spectrum_spin_rate) = self.stellar_spectrum_spin_rate {
            imaginary_k2 *= (spin_rate / spectrum_spin_rate).powi(2);
        }

        (real_k2, imaginary_k2)
    }

    // Find the real part and the imaginary part of the Love number associated with
    // the excitation frequenccy wk2
    fn interpolate_love_numbers(&self, wk2: f64) -> (f64, f64) {
        let im_k2;
        let re_k2;
        // Find the index of the closest match to use for the interpolation
        if wk2 <= self.spectrum_excitation_frequency[0] {
            // If wk2 is less than or equal to the first element, take the first value
            im_k2 = self.spectrum_imaginary_part[0];
            re_k2 = self.spectrum_real_part[0];
        } else if wk2
            >= self.spectrum_excitation_frequency[self.spectrum_excitation_frequency.len() - 1]
        {
            // If wk2 is greater than or equal to the last element, take the last value
            im_k2 = self.spectrum_imaginary_part[self.spectrum_excitation_frequency.len() - 1];
            re_k2 = self.spectrum_real_part[self.spectrum_excitation_frequency.len() - 1];
        } else {
            // Find the index of the closest match to use for the interpolation
            match self
                .spectrum_excitation_frequency
                .binary_search_by(|val| val.total_cmp(&wk2))
            {
                Ok(i) => {
                    // Exact match found: love_number[i] == wk2
                    assert!(i < self.spectrum_real_part.len());
                    re_k2 = self.spectrum_real_part[i];
                    assert!(i < self.spectrum_imaginary_part.len());
                    im_k2 = self.spectrum_imaginary_part[i];
                }
                Err(i) => {
                    // wk2 is between love_number[i - 1] and love_number[i]
                    assert!(i < self.spectrum_real_part.len());
                    assert!(i < self.spectrum_imaginary_part.len());
                    let prev_freq = self.spectrum_excitation_frequency[i - 1];
                    let next_freq = self.spectrum_excitation_frequency[i];
                    let delta = (wk2 - prev_freq) / (next_freq - prev_freq);
                    im_k2 = (1.0 - delta) * self.spectrum_imaginary_part[i - 1]
                        + delta * self.spectrum_imaginary_part[i];
                    re_k2 = (1.0 - delta) * self.spectrum_real_part[i - 1]
                        + delta * self.spectrum_real_part[i];
                }
            }
        }

        (re_k2, im_k2)
    }
}
