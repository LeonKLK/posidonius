# Known Bugs: CTL Tidal Model + Stellar/Planetary Evolution

Two pre-existing bugs found in the coupling between the Constant Time Lag (CTL)
tidal model and the evolution models. Both verified against commit `491c643`
("cleanup, refactor, optimise (#3)").

Status: **documented, not yet fixed** (2026-08-08).

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
