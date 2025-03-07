use serde::{Deserialize, Serialize};

#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct Polynomials {
    pub eccentricity_function_g_2pq: [[f64; 15]; 3],
    pub eccentricity_function_g_3pq: [[f64; 15]; 4],
    pub inclination_function_f_2mp: [[f64; 3]; 3],
    pub inclination_function_f_3mp: [[f64; 4]; 4],
}

impl Polynomials {
    // Constructor function that creates a Polynomials instance with all values set to 0.0
    pub fn new() -> Self {
        Polynomials {
            eccentricity_function_g_2pq: [[0.0; 15]; 3],
            eccentricity_function_g_3pq: [[0.0; 15]; 4],

            inclination_function_f_2mp: [[0.0; 3]; 3],
            inclination_function_f_3mp: [[0.0; 4]; 4],
        }
    }
    pub fn update_inclination(&mut self, obliquity: f64) {
        self.calculate_inclination_function_f_2mp(obliquity);
        self.calculate_inclination_function_f_3mp(obliquity);
    }

    pub fn update_eccentricity_2d(&mut self, eccentricity: f64) {
        self.calculate_eccentricity_function_g_2pq(eccentricity);
    }

    pub fn update_eccentricity_3d(&mut self, eccentricity: f64) {
        self.calculate_eccentricity_function_g_3pq(eccentricity);
    }

    // The eccentricity function of the Kaula 1961 developpement Glpq(e) from the Hensen coeff [?] (see Cayley table [?])
    fn calculate_eccentricity_function_g_2pq(&mut self, eccentricity: f64) {
        // in this notation: [0, 1, 2, ..., 7, ..., 12, 13, 14] maps to [-7, -6, -5, ..., 0, ..., 5, 6, 7]
        let ecc = eccentricity;
        let ecc_2 = ecc * ecc;
        let ecc_3 = ecc_2 * ecc;
        let ecc_4 = ecc_3 * ecc;
        let ecc_5 = ecc_4 * ecc;
        let ecc_6 = ecc_5 * ecc;
        let ecc_7 = ecc_6 * ecc;

        let borrowed_slice = &mut self.eccentricity_function_g_2pq[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = (15_625. / 129_024.) * ecc_7;
        first[1] = (4. / 45.) * ecc_6;
        first[2] = (81. / 1280.) * ecc_5 + (81. / 2048.) * ecc_7;
        first[3] = (1. / 24.) * ecc_4 + (7. / 240.) * ecc_6;
        first[4] = (1. / 48.) * ecc_3 + (11. / 768.) * ecc_5 + (313. / 30_720.) * ecc_7;
        first[5] = 0.;
        first[6] = -0.5 * ecc + (1. / 16.) * ecc_3 - (5. / 384.) * ecc_5 - (143. / 18_432.) * ecc_7;
        first[7] = 1. - (5. / 2.) * ecc_2 + (13. / 16.) * ecc_4 - (35. / 288.) * ecc_6;
        first[8] = (7. / 2.) * ecc - (123. / 16.) * ecc_3 + (489. / 128.) * ecc_5
            - (1763. / 2048.) * ecc_7;
        first[9] = (17. / 2.) * ecc_2 - (115. / 16.) * ecc_4 + (601. / 48.) * ecc_6;
        first[10] = (845. / 48.) * ecc_3 - (32_525. / 768.) * ecc_5 + (208_225. / 6144.) * ecc_7;
        first[11] = (533. / 16.) * ecc_4 - (13_827. / 160.) * ecc_6;
        first[12] = (228_347. / 3840.) * ecc_5 - (3_071_075. / 18_432.) * ecc_7;
        first[13] = (73_369. / 720.) * ecc_6;
        first[14] = (12_144_273. / 71_680.) * ecc_7;

        second[0] = (432_091. / 30_720.) * ecc_7;
        second[1] = (3167. / 320.) * ecc_6;
        second[2] = (1773. / 256.) * ecc_5 - (4987. / 6144.) * ecc_7;
        second[3] = (77. / 16.) * ecc_4 + (129. / 160.) * ecc_6;
        second[4] = (53. / 16.) * ecc_3 + (393. / 256.) * ecc_5 + (24_753. / 10_240.) * ecc_7;
        second[5] = (9. / 4.) * ecc_2 + (7. / 4.) * ecc_4 + (141. / 64.) * ecc_6;
        second[6] = (3. / 2.) * ecc
            + (27. / 16.) * ecc_3
            + (261. / 128.) * ecc_5
            + (14_309. / 6144.) * ecc_7;
        second[7] = (1. - ecc_2).powf(-3. / 2.);
        // array reflected here
        let (left, right) = second.split_at_mut(8);
        right.copy_from_slice(&left[0..=6]);
        right.reverse();

        // Copy in reverse order
        third.copy_from_slice(&first[..]);
        third.reverse();
    }

    // The eccentricity function of the Kaula 1961 developpement Glpq(e) from the Hensen coeff [?] (see Cayley table [?])
    fn calculate_eccentricity_function_g_3pq(&mut self, eccentricity: f64) {
        // in this notation: [0, 1, 2, ..., 7, ..., 12, 13, 14] maps to [-7, -6, -5, ..., 0, ..., 5, 6, 7]
        let ecc = eccentricity;
        let ecc_2 = ecc.powi(2);
        let ecc_3 = ecc.powi(3);
        let ecc_4 = ecc.powi(4);
        let ecc_5 = ecc.powi(5);
        let ecc_6 = ecc.powi(6);
        let ecc_7 = ecc.powi(7);

        let borrowed_slice = &mut self.eccentricity_function_g_3pq[..];
        let [first, second, third, fourth] = borrowed_slice else {
            unreachable!()
        };

        first[0] = (8. / 315.) * ecc_7;
        first[1] = (81. / 5120.) * ecc_6;
        first[2] = (1. / 120.) * ecc_5 + (13. / 1440.) * ecc_7;
        first[3] = (1. / 384.) * ecc_4 + (1. / 384.) * ecc_6;
        first[4] = 0.;
        first[5] = (1. / 8.) * ecc_2 + (1. / 48.) * ecc_4 + (55. / 3072.) * ecc_6;
        first[6] = -ecc + (5. / 4.) * ecc_3 - (7. / 48.) * ecc_5 + (23. / 288.) * ecc_7;
        first[7] = 1. - 6. * ecc_2 + (423. / 64.) * ecc_4 - (125. / 64.) * ecc_6;
        first[8] = 5. * ecc - 22. * ecc_3 + (607. / 24.) * ecc_5 - (98. / 9.) * ecc_7;
        first[9] = (127. / 8.) * ecc_2 - (3065. / 48.) * ecc_4 + (243_805. / 3072.) * ecc_6;
        first[10] = (163. / 4.) * ecc_3 - (2577. / 16.) * ecc_5 + (1089. / 5.) * ecc_7;
        first[11] = (35_413. / 384.) * ecc_4 - (709_471. / 1920.) * ecc_6;
        first[12] = (23_029. / 120.) * ecc_5 - (35_614. / 45.) * ecc_7;
        first[13] = (385_095. / 1024.) * ecc_6;
        first[14] = (44_377. / 63.) * ecc_7;

        second[0] = (16_337. / 2240.) * ecc_7;
        second[1] = (48_203. / 9240.) * ecc_6;
        second[2] = (899. / 240.) * ecc_5 + (2441. / 480.) * ecc_7;
        second[3] = (343. / 128.) * ecc_4 + (2819. / 640.) * ecc_6;
        second[4] = (23. / 12.) * ecc_3 + (89. / 24.) * ecc_5 + (5663. / 960.) * ecc_7;
        second[5] = (11. / 8.) * ecc_2 + (49. / 16.) * ecc_4 + (15_665. / 3072.) * ecc_6;
        second[6] = ecc * (1. - ecc_2).powf(-5. / 2.);
        second[7] = 1. + (1. / 2.) * ecc_2 + (239. / 64.) * ecc_4 - (3323. / 576.) * ecc_6;
        second[8] = 3. * ecc + (11. / 4.) * ecc_3 + (245. / 48.) * ecc_5 + (463. / 64.) * ecc_7;
        second[9] = (53. / 8.) * ecc_2 + (39. / 16.) * ecc_4 + (7041. / 1024.) * ecc_6;
        second[10] = (163. / 4.) * ecc_3 - (2577. / 16.) * ecc_5 + (1089. / 5.) * ecc_7;
        second[11] = (35_413. / 384.) * ecc_4 - (709_471. / 1920.) * ecc_6;
        second[12] = (23_029. / 120.) * ecc_5 - (35_614. / 45.) * ecc_7;
        second[13] = (385_095. / 1024.) * ecc_6;
        second[14] = (44_377. / 63.) * ecc_7;

        // Copy in reverse order
        third.copy_from_slice(&second[..]);
        third.reverse();

        // Copy in reverse order
        fourth.copy_from_slice(&first[..]);
        fourth.reverse();
    }

    // The inclination function of the Kaula 1961 developpement Flmp(i)
    fn calculate_inclination_function_f_2mp(&mut self, inclination: f64) {
        let sin_inc = sin!(inclination);
        let cos_inc = cos!(inclination);
        let sin_inc_2 = sin_inc * sin_inc;

        let borrowed_slice = &mut self.inclination_function_f_2mp[..];
        let [first, second, third] = borrowed_slice else {
            unreachable!()
        };

        first[0] = -(3. / 8.) * sin_inc_2;
        first[1] = (3. / 4.) * sin_inc_2 - (1. / 2.);
        first[2] = -(3. / 8.) * sin_inc_2;

        second[0] = (3. / 4.) * sin_inc * (1. + cos_inc);
        second[1] = -(3. / 2.) * sin_inc * cos_inc;
        second[2] = (3. / 4.) * sin_inc * (cos_inc - 1.);

        third[0] = (3. / 4.) * (1. + cos_inc).powi(2);
        third[1] = (3. / 2.) * sin_inc_2;
        third[2] = (3. / 4.) * (1. - cos_inc).powi(2);
    }

    // F3pmp
    // The inclination function of the Kaula 1961 developpement Flmp(i)
    pub fn calculate_inclination_function_f_3mp(&mut self, inclination: f64) {
        let sin_inc = f64::sin(inclination);
        let cos_inc = f64::cos(inclination);
        let cos_inc_2 = cos_inc.powi(2);
        let sin_inc_2 = sin_inc.powi(2);
        let sin_inc_3 = sin_inc.powi(3);

        let borrowed_slice = &mut self.inclination_function_f_3mp[..];
        let [first, second, third, fourth] = borrowed_slice else {
            unreachable!()
        };

        let fac_a = (5. / 16.) * sin_inc_3;
        let fac_b = (3. / 4.) * sin_inc;

        first[0] = -fac_a;
        first[1] = fac_a - fac_b;
        first[2] = fac_b - fac_a;
        first[3] = fac_a;

        let fac_c = (15. / 16.) * sin_inc_2;
        let fac_d = fac_c - (3. / 4.);
        let fac_e = 3.0 * cos_inc;

        second[0] = -fac_c * (1.0 + fac_e);
        second[1] = fac_d * (1.0 + fac_e);
        second[2] = fac_d * (1.0 - fac_e);
        second[3] = -fac_c * (1.0 - fac_e);

        let fac_f = (15. / 8.) * sin_inc_2;
        let fac_h = 2.0 * cos_inc;
        let fac_k = 3.0 * cos_inc_2;

        third[0] = fac_f * (1.0 + cos_inc).powi(2);
        third[1] = fac_f * (1.0 - fac_h - fac_k);
        third[2] = -fac_f * (1.0 + fac_h - fac_k);
        third[3] = -fac_f * (1.0 - cos_inc).powi(2);

        let fac_x = 15. / 8.;
        let fac_y = (48. / 8.) * sin_inc_2;

        fourth[0] = fac_x * (1.0 + cos_inc).powi(3);
        fourth[1] = fac_y * (1.0 + cos_inc);
        fourth[2] = fac_y * (1.0 - cos_inc);
        fourth[3] = fac_x * (1.0 - cos_inc).powi(3);
    }
}
