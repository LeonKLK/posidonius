# Known Bugs: Tidal Models (CTL + Kaula) and Evolution

Pre-existing bugs found in the tidal-effect code. All verified against commit
`491c643` ("cleanup, refactor, optimise (#3)").

Status (2026-08-10):

| # | Bug | Affected configuration | Status |
|---|-----|------------------------|--------|
| 1 | Love-number evolution is a silent no-op | LeconteChabrier2013 + CTL on the same body | **Open** (low priority, unused config) |
| 2 | Asymmetric HashMap key in dynamical-tide override removal | Orbiting body with BM2016/GB2017/LC2013(true) + CTL | **Fixed** in the tidal restructure (remove key now matches the set key) |
| 3 | Kaula stellar torque uses only the last planet's force | Kaula star (CentralBody) + 2 or more planets | **Fixed** in the tidal restructure (stellar secular force stored per companion; order-invariance covered by a regression test) |

---

## Bug 1 — Love-number evolution is a silent no-op

### Location
- `src/effects/evolution.rs:549-557` (function `calculate_particles_non_spin_dependent_evolving_quantities`)

### Affected configuration
- **Only observable with:** `EvolutionType::LeconteChabrier2013` + `TidalModel::ConstantTimeLag`
  on the same body (evolving gas giant / brown dwarf with CTL tide).
- Why only this combination:
  1. `Evolver::love_number()` (`evolution.rs:473-484`) returns a *new* interpolated
     value only for `LeconteChabrier2013`; for every other evolution type it returns
     the input unchanged, so a broken write-back has no observable effect.
  2. The write-back is gated on `TidalModel::ConstantTimeLag`. Kaula uses its
     frequency spectrum and CreepCoplanar its own parameters — the scalar
     `love_number` only feeds the CTL radial force.
- Note: the buggy *pattern* would silently swallow any future evolution model that
  evolves k2, so it is worth fixing regardless.

### Root cause
`TidalModel` and `ConstantTimeLagParameters` both derive `Copy`
(`src/effects/tides.rs:89-95`, `src/effects/tides/constant_time_lag.rs:10`).
The match pattern binds **by copy**, not by reference:

```rust
// evolution.rs:549-557
match &mut particle.tides.effect {
    &mut (TidesEffect::CentralBody(mut tidal_model)
    | TidesEffect::OrbitingBody(mut tidal_model)) => {   // <- copies out of the &mut
        if let &mut TidalModel::ConstantTimeLag(mut params) = &mut tidal_model {
            // ^ copies AGAIN out of the local copy
            params.love_number = evolver.love_number(current_time, params.love_number);
            // ^ mutates a stack-local temporary, dropped at end of block
        }
    }
    TidesEffect::Disabled => {}
}
```

The pattern `&mut (Variant(mut binding))` matches *through* the mutable reference
and, because the inner type is `Copy`, `tidal_model` is a stack-local duplicate.
The `if let &mut ...(mut params)` copies a second time. The assignment mutates
`params` (a copy of a copy); nothing is ever written back to
`particle.tides.effect`.

Contrast with the radius/rg² update just above (`evolution.rs:534-544`), which
correctly assigns through the real reference (`particle.radius = new_radius;`).

### Consequence
For an evolving LeconteChabrier2013 body with a CTL tide:
- The k2 column of the evolution table is loaded and linearly interpolated every
  timestep — then thrown away.
- The CTL radial force (`constant_time_lag.rs:439-450`) reads `params.love_number`
  from `particle.tides.effect`, which forever holds the **initial** value.
- The Python docs (`posidonius/effects/evolution.py:32-35`) tell users the initial
  love number "will be ignored" under Jupiter evolution — in reality the initial
  value is the *only* one ever used.

### Suggested fix
Bind mutably by reference instead of by copy:

```rust
match &mut particle.tides.effect {
    TidesEffect::CentralBody(tidal_model)
    | TidesEffect::OrbitingBody(tidal_model) => {
        if let TidalModel::ConstantTimeLag(params) = tidal_model {
            params.love_number = evolver.love_number(current_time, params.love_number);
        }
    }
    TidesEffect::Disabled => {}
}
```

(No `&mut` in the patterns, no `mut` on the bindings — match ergonomics then give
`tidal_model: &mut TidalModel` and `params: &mut ConstantTimeLagParameters`, so
the assignment writes through to the particle.)

### Verification idea
Unit/integration test: build a body with `LeconteChabrier2013` evolution + CTL
tide, step the evolver past a table point where k2 changes, then assert that the
love number stored in `particle.tides.effect` differs from the initial value.
Currently such a test would fail.

---

## Bug 2 — Asymmetric HashMap key in dynamical-tide override removal

### Location
- `src/effects/tides/constant_time_lag.rs:185-197`
  (function `calculate_host_dependent_scaled_dissipation_factors`)
- Correct sibling for comparison: `constant_time_lag.rs:110-122`
  (function `calculate_planet_dependent_scaled_dissipation_factors`)

### Affected configuration
- **Only fires when an ORBITING body carries** `BolmontMathis2016`,
  `GalletBolmont2017`, or `LeconteChabrier2013(true)` evolution **with a CTL tide**
  (gate at `constant_time_lag.rs:139-145`).
- BM2016 / GB2017 are stellar tables, so in practice this means an evolving
  gas giant with `LeconteChabrier2013(true)` — the "evolving hot Jupiter with
  dynamical tides" scenario.
- **Latent (harmless) for the usual setup** where only the star (tidal host)
  evolves with these models: then only the correct sibling function runs.

### Background
When the frequency-averaged dynamical tide is active, per-pair dissipation
overrides live in `Universe.pair_dependent_scaled_dissipation_factor`
(`src/particles/universe.rs:61`), keyed directionally:

```
key = id * MAX_PARTICLES + depends_on_id
// "dissipation of body `id`, as modified by companion `depends_on_id`"
```

Every force evaluation, two symmetric functions refresh the map
(`constant_time_lag.rs:18-36`):
1. `calculate_planet_dependent_...` — the STAR's dissipation per planet,
   entry key `(star, planet)`.
2. `calculate_host_dependent_...` — the evolving PLANET's own dissipation,
   entry key `(planet, star)`.

Each has an if/else on the corotation gate `|Ω − n| < Ω`:
dynamical regime → `set` the override; equilibrium regime → `remove` it so the
force falls back to the static equilibrium sigma
(`get_pair_dependent_scaled_dissipation_factor_or_else`, lines 226-245).

### Root cause
In function 2, `set` and `remove` use **different keys**:

```rust
// constant_time_lag.rs — set (dynamical regime): key = (particle, host)  ✓
set_pair_dependent_scaled_dissipation_factor(
    pair_dependent_scaled_dissipation_factor,
    particle.id,
    tidal_host_particle.id,
    host_dependent_scaled_dissipation_factor,
);
} else {
// remove (equilibrium regime): key = (host, particle)  ✗ — the STAR's entry
remove_pair_dependent_scaled_dissipation_factor(
    pair_dependent_scaled_dissipation_factor,
    tidal_host_particle.id,
    particle.id,
);
```

The sibling function 1 is consistent (both `set` and `remove` use
`(star, planet)`, lines 110-122).

### Consequence
Both failure modes fire on the timestep the evolving planet crosses from the
dynamical into the equilibrium regime (`|Ω_p − n| ≥ Ω_p`):

1. **Stale value never cleared.** The planet's entry `(planet, star)` stays in
   the map. Readers keyed `(planet, star)` keep finding the frozen dynamical
   dissipation from the last dynamical-regime step — the planet dissipates at a
   stale (typically much larger) rate indefinitely instead of falling back to
   its equilibrium sigma.
2. **The star's entry is wrongly deleted.** The removed key `(star, planet)` is
   exactly what function 1 just set (function 1 runs first, lines 24-35). A star
   legitimately in the dynamical regime silently loses its override and drops
   back to equilibrium dissipation — invisible in the output.

### Suggested fix
Swap the two arguments in the `remove` call so it matches the `set`:

```rust
remove_pair_dependent_scaled_dissipation_factor(
    pair_dependent_scaled_dissipation_factor,
    particle.id,
    tidal_host_particle.id,
);
```

### Verification idea
Test with a two-body system where the orbiting body has
`LeconteChabrier2013(true)` + CTL and spin/orbit chosen so it starts in the
dynamical regime, then transitions out (e.g. spin down the planet). Assert
after the transition that:
- the map no longer contains key `planet.id * MAX_PARTICLES + star.id`, and
- the star's key `star.id * MAX_PARTICLES + planet.id` is still present when the
  star remains in the dynamical regime.

---

## Bug 3 — Kaula stellar torque uses only the last planet's force (multi-planet)

### Location
- Storage slot (single, on the star): `src/effects/tides/kaula.rs:141`
  (`particle.tides.get_kaula_mut().tidal_force = secular_projection;` inside
  `calculate_tidal_force_component`)
- Overwriting loop: `src/effects/tides.rs:503-529` (stellar-tide loop inside
  `calculate_tidal_acceleration`, gated on the star being
  `CentralBody(TidalModel::Kaula)`)
- Corrupted consumer: `src/effects/tides.rs:347-365` (star branch of
  `calculate_dangular_momentum_dt_due_to_tides`) →
  `src/effects/tides/kaula.rs:253-288` (`calculate_torque_due_to_tides`,
  `central_body == true` branch reads the star's single stored force at line 269)

### Affected configuration
- **Kaula stellar tide (`CentralBody(Kaula)`) with 2 or more planets.**
- Single-planet systems are **exact** — the sum has one term and the stored
  force belongs to that planet.
- Independent of what tidal model the planets carry (CTL, Kaula, or none).

### Root cause
`KaulaParameters` has a single `tidal_force` field. The stellar acceleration
loop iterates the planets in array order and, for each planet `i`, computes the
correct stellar secular force `F_i` — but stores it in that one slot on the
star, overwriting the previous iteration. After the loop the slot holds only
`F_N`, the force due to the **last planet in the particle array**.

The star-torque loop then iterates planets again and, for each planet, calls
the Kaula torque function, which reads the star's single stored force:

```text
dL/dt (star) = - Σ_i  r_i × F_stored  =  - Σ_i  r_i × F_N     (WRONG)
                                          correct: Σ_i r_i × F_i
```

Each planet contributes its own position vector `r_i`, but every cross product
uses the last planet's force `F_N`. Only the `i = N` term is physically
correct; the other `N-1` terms mix one planet's geometry with another planet's
force amplitude.

### Consequence
- **Planetary orbits remain correct** — inside the stellar acceleration loop
  each iteration applies its own `F_i` before the slot is overwritten.
- **The star's spin evolution (`dangular_momentum_dt`) is wrong** in both
  magnitude and potentially direction. Because the tidal force scales steeply
  with distance, if the last planet in the array is close-in, the outer
  planets' pseudo-torques are inflated by the inner planet's force amplitude.
- The error **feeds back into the forces over time**: the stellar spin enters
  the tidal excitation frequencies `w_2mpq = (2-2p+q)·n - m·Ω` used in the
  love-number interpolation of subsequent timesteps.
- Note: the same storage pattern on a Kaula *planet* is harmless — each planet
  owns its own slot, written and read one-to-one.

### Suggested fix
Store the stellar secular force **per planet** instead of once on the star.
Two natural options:
1. Add a per-planet field in `particle.tides.parameters.internal` (mirroring
   how CTL stores `orthogonal_component_of_the_tidal_force_due_to_stellar_tide`
   per planet), and have the star-torque branch read the field from the
   *planet* being iterated; or
2. Compute the torque contribution directly inside the stellar acceleration
   loop (`tides.rs:508-528`) where `F_i` is still in hand, accumulating
   `Σ r_i × F_i` there, and skip the Kaula star branch in the torque loop.

Option 2 fits naturally into the planned restructuring of the stellar-tide
loop (decoupling stellar tide from the planet's tidal model).

### Verification idea
Two-planet system with a Kaula star: compute the stellar
`dangular_momentum_dt` once with planets in order [A, B] and once with
[B, A]. Currently the results differ (the stored force belongs to whichever
planet comes last); after the fix they must be identical and equal to
`-(r_A × F_A + r_B × F_B)`.

---

## Related context (not bugs, design notes)

- The dynamical tide is *added* to the equilibrium `dissipation_factor`
  (`constant_time_lag.rs:103-108`), and its activation is implicit in the choice
  of evolution table — there is no independent switch in the CTL configuration.
  A planned restructuring will separate "CTL equilibrium tide" from "CTL
  dynamical tide".
- A commented-out `panic!` at `constant_time_lag.rs:126` and `:204` warns these
  evolution models "may not be ready yet for scientific exploitation".
- There is no consistency check (Rust or Python) that a body using
  BM2016/GB2017/LC2013(true) actually uses the CTL model; with Kaula or
  CreepCoplanar the dynamical-tide machinery silently does nothing.
