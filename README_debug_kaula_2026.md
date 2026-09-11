# debug_kaula_2026 — Kaula tidal model: known problems and inconsistencies

Branch created from `master` (= DynaClim `491c643`, "cleanup, refactor, optimise (#3)") on 2026-09-11.
Purpose: make the Kaula stellar tide in Posidonius consistent with the secular code Spiroid
(1 star + 1 planet, star tides only) and fix what is found on the way.

Evidence for the findings below comes from runs stored in `~/Documents/spi_pos_comparison`
(case `cases/a03_k2_flatline_signed.py`, 1 Msun + 5 Mearth at 0.03 AU, flat signed Im(k2) spectrum
shared with Spiroid, WHFast, dt = 0.01 d) and from the analytic rate
`T = (3/2) Im(k2) G m_p^2 R_*^5 / a^6` for a coplanar circular orbit.

## Status legend
- [ ] unresolved  - [x] fixed on this branch (see commit)

## Findings

### F1. [ ] Exactly circular orbit gives ZERO Kaula tidal force (star and planet tides alike)
`src/effects/tides/kaula/components_2d.rs::zero_eccentricity_components` (taken when
`keplerian_elements.eccentricity == 0.0`) reads `love_numbers.real(0,1,0)` / `imaginary(2,0,0)`.
The last argument is the cache INDEX q (0..14, index 7 = physical q = 0), so index 0 means physical
q = -7. `LoveNumber::refresh_cache` only fills the indices returned by `select_eccentricty_order_q`
(5..9 for e <= 0.1 in `491c643`), so index 0 is never written below e > 0.30 and stays 0.
Evidence: tides ON and OFF are bit-identical for e = 0 (3 d and 100 d). With e = 1e-6 the tide acts.
The repository tests never see it because every test case uses e = 1e-6 or e = 0.1.
Fix direction: read index 7 (physical q = 0), or route e == 0 through `calculate_2d_components`.
Related: the local branch `zero_ecc_in_kaula` (June 2026) contains an earlier attempt.

### F2. [x] Stellar (central-body) tidal force built with the planet-tide prefactor
Fixed by the role-based refactor (commit "kaula: name inputs by role"); references regenerated in the
following commit. Verified: da/dt = -6.906e-7 m/s on the comparison case, as the analytic rate.
For the stellar tide, `tides.rs::calculate_tidal_acceleration` calls
`kaula::calculate_tidal_force(planet, star, central_body = true)`, so inside `kaula.rs`
`tidal_host_particle` is the PLANET and `particle` is the STAR.
`calculate_base_constant(first, second, r, a) = G * second.mass^2 * first.radius^5 / (a^6 r)`,
i.e. the first argument must be the tidally deformed body.
- `components_2d.rs::orthogonal_constant`, central branch: passes (planet, star) -> `G M_*^2 R_p^5`. WRONG.
- `components_2d.rs::radial_constant`, central branch: passes (star, planet) -> `G m_p^2 R_*^5`. Correct.
- `components_3d.rs` normal (l.345), orthogonal (l.599) and radial (l.635) central branches: all pass
  (planet, star). WRONG in every 3D component.
Consequence: in a coplanar circular orbit the tangential force sets BOTH the orbital migration and the
stellar spin torque, so both are off by `(M_*/m_p)^2 (R_p/R_*)^5` (0.37 for the 5 Mearth case:
measured da/dt = -2.5506e-7 m/s vs -6.906e-7 expected; dOmega/dt = 4.361e-24 vs 1.181e-23 rad/s^2;
doubling R_p with the planet tide disabled multiplies the stellar-tide rates by exactly 32 = 2^5).
Planetary tides (non-central branches) use the correct order everywhere.
Fix direction: give the constants explicit roles (`tidal_deformed_body`, `tidal_perturber`, orbit
quantities always from the planet). Note: swapping the two arguments of `orthogonal_constant` alone
would give NaN, because `calculate_2d_constant` takes sin(theta) from the FIRST argument's
`tides.coordinates.position`, which is zeroed for the host in `inertial_to_heliocentric_coordinates`.
After the fix the stored reference data of `enabled_star_tides` and `enabled_both_tides`
(`tests/data/test_tides_kaula-*`) must be regenerated: they were produced with the wrong prefactor.

### F3. [ ] 2D vs 3D sign convention of the heliocentric distance for the central body
(After the refactor the 2D sign flip lives in `components_2d.rs::signed_heliocentric_distance`; the 3D
constants take the unsigned `orbit.heliocentric_distance` for both tides, as before.)
2D central branches pass `-tidal_host_particle.heliocentric_distance` to the constants, 3D central
branches pass `+tidal_host_particle.heliocentric_distance` (normal, radial) or `1.0` (orthogonal).
Not yet checked whether 3D compensates the sign elsewhere. The 3D stellar tide has no test.

### F4. [ ] Circular-branch radial force: `imk2_2200` vs `rek2_2200`
Uncommitted comment found in the working tree (kept on this branch):
`// Bug exist for now: imk2_2200 should be rek2_2200` above `zero_eccentricity_components`.
The 9/4 term of the radial force uses Im(k2_2200); the general branch `calculate_2d_components` uses
`rek2_220q` for the same term (`sum_over_j_3`), so the circular branch is inconsistent with it.
Only matters once F1 is fixed (today the branch returns 0 anyway).

### F5. [ ] `select_eccentricty_order_q` edit breaks the planetary-tide reference test
An uncommitted edit (saved as `~/Documents/spi_pos_comparison/posidonius_uncommitted_q_range_edit.patch`,
NOT applied on this branch) replaced the `e > 0.10 -> (4,11)` / `e <= 0.1 -> (5,10)` tiers by
`e > 1e-8 -> (5,10)`, `e > 0 -> (6,9)`, `e == 0 -> (7,8)`. The test planet starts at e = 0.1, i.e. on
the boundary of the removed 0.10..0.15 tier, so `enabled_planet_tides_rust` and `enabled_both_tides_rust`
fail (positions differ from the stored reference at 1e-14). Any change of the q ranges is a physics
change and needs new reference data, not a silent edit.

### F6. [x] Missing test input file (local only, `input/` is gitignored)
`tests/test_tides_kaula.rs` needs
`input/love_numbers/Results_Trappist1_h_Fe_90_Si_2_170K_visco_freq_Imk2_posidonius_for_testing.txt`.
The whole `input/` directory is in `.gitignore` (data is distributed separately), so the file cannot be
committed; without it 4 of 7 Kaula tests panic with NotFound. A copy was placed locally from
`spiroid/comparison/alex_venus/` (471 rows, identical to the arrays stored in
`tests/data/test_tides_kaula-enabled_planet_tides/case.json`). With it, clean `491c643` passes all 7.
Keep that file in `input/love_numbers/` on every machine that runs the tests.

### F7. [ ] Kaula tests overflow the default test-thread stack
`cargo test --test test_tides_kaula` aborts with "stack overflow" when run multi-threaded
(`LoveNumber` holds 3 x 1024 f64 by value, `Particle` copies are large). Run with
`RUST_MIN_STACK=268435456 cargo test --release --test test_tides_kaula -- --test-threads=1`.

## Other differences vs Spiroid worth remembering (not bugs)
- Stellar evolution: Posidonius `GalletBolmont2017` gives R = 1.4917 Rsun at 5 Myr; Spiroid Starevol
  table (`savgol_10.csv`) gives 1.4445 Rsun. Tidal torque scales as R^5 (17 % difference).
- Spiroid is a two-zone star (tidal torque on the convective envelope, core-envelope coupling);
  Posidonius spins the star as one body with `radius_of_gyration_2`.
- Wind prescriptions differ (Bouvier 1997 here vs Matt+2015 in Spiroid).
- The comparison case sits ~0.5 % from corotation (2n - 2Omega = 0.009 spin_spec), so the SIGN of the
  stellar torque depends on tiny differences in Omega(t).

## Role convention after the refactor (kaula.rs)
- `tidal_deformed_body`: star for the stellar tide (`central_body == true`), planet for the planetary
  tide. Provides spin, radius and the love number spectrum.
- `tidal_perturber`: the other body. Provides only its mass.
- `Orbit::from_planet(planet)`: heliocentric position/velocity/distance and the `tides.coordinates`
  position/distance of the PLANET, used for the keplerian elements and the projection angles whichever
  body is deformed.
- `central_body` only selects: love-number parity, the slot where the orthogonal component is stored
  on the planet, the sign of the heliocentric distance in the 2D prefactors, the sign of r in the torque.

## Shared stellar profile with Spiroid: `EvolutionType::Starevol`
`scripts/make_starevol_profile.py` converts Spiroid's `examples/data/star/evolution/savgol_10.csv` into
`input/Starevol/M_10.dat` (age[yr], radius[Rsun], rg2 = (I_rad + I_conv)/(M R^2), I_conv/(M R^2), mass). Rust
`EvolutionType::Starevol(mass)` and python `posidonius.Starevol(mass)` read it (1 Msun only, mass accepted in
0.95..1.05 so the star can be given Spiroid's SOLAR_MASS in kg). `input/` is gitignored: regenerate the file locally.
The pip-installed package must be reinstalled (`pip install .`) for `cases/*.py` to see new python classes.
Result (comparison case, e = 1e-6, 100 d): with the same profile Posidonius and Spiroid give da/dt = -5.7477e-7 and
-5.7478e-7 m/s (ratio 1.00002), both equal to the analytic rate. With GalletBolmont2017 the ratio was 0.832 = (R ratio)^5.
The stellar spin still differs: Spiroid applies the torque to the convective envelope (two-zone star) and has a
different wind; Posidonius spins the whole star with rg2_total. Not a tide issue.

## How to run the reference tests
```
RUST_MIN_STACK=268435456 cargo test --release --test test_tides_kaula -- --test-threads=1
```
`enabled_planet_tides_rust` is the guard for refactors: the planetary Kaula tide is believed correct
and its stored positions must stay bit-identical (1e-14).
