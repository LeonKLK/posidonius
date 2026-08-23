# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Posidonius is an N-body simulator for planetary systems with tidal effects, rotational flattening, general relativity and stellar/planetary evolution. It is a dual-language project:

- **Rust crate** (`src/`): the actual simulator, installed as the `posidonius` binary (`posidonius start case.json snapshot.bin history.bin`, `posidonius resume ...`).
- **Python package** (`posidonius/`): generates simulation cases as JSON files that the Rust binary consumes, and analyzes the binary history output.

The two sides communicate through a **JSON contract**: the Python classes in `posidonius/effects/*.py` and `posidonius/particles/*.py` build nested dicts that must deserialize field-for-field into the Rust serde structs (`src/effects/*.rs`, `src/particles/*.rs`). **Any change to a Rust struct that is serialized must be mirrored in the corresponding Python dict** (and vice versa). New Rust fields should carry `#[serde(default)]` so older JSON case files keep loading. Recovery snapshots (`.bin`) are bincode and are NOT tolerant to struct layout changes.

## Commands

```bash
# Build / install the simulator
cargo build --release
cargo install --path . --force

# Rust tests — RUST_MIN_STACK is required (Kaula uses large stack arrays)
RUST_MIN_STACK=33554432 cargo test
RUST_MIN_STACK=33554432 cargo test --lib            # fast: unit + golden regression tests only
RUST_MIN_STACK=33554432 cargo test --test test_tides enabled_tides_rust   # single test
cargo fmt && cargo clippy

# Python tests (needs `pip install '.[dev]'`)
pytest

# Benchmarks (criterion): baseline before a change, compare after
cargo bench -- --save-baseline reference
cargo bench -- --baseline reference
```

### External input data

Most integration tests and the Python package require evolution tables and love-number spectra that are NOT in the repo:

```bash
curl -O https://www.blancocuaresma.com/s/repository/posidonius/input.tar.gz
tar -zxvf input.tar.gz && rm -f input.tar.gz   # creates ./input/
```

Without `input/`: `import posidonius` raises at import time, and most tests in `tests/` fail in `Evolver::new` (file not found). The **golden regression tests** (`cargo test --lib golden`, in `src/effects/tides/golden_tests.rs`) are deliberately self-contained (synthetic spectra, non-evolving bodies) and always run — they freeze the numerical output of the tidal pipeline and are the primary safety net for refactors of the tides code.

### Reference-data regeneration

Rust integration tests compare against stored JSON references in `tests/data/`. If a numerical difference is *intended*:

```bash
find tests/data/ -name 'particle_*.json' -delete && RUST_MIN_STACK=33554432 cargo test
```

Python tests similarly regenerate via `find posidonius/tests/data/ -name 'case.json' -delete && pytest`. Rust-generated and Python-generated `case.json` files are compared by the `*_rust_vs_python` tests; `scripts/clean_json.py` reformats the Rust ones for comparability.

## Architecture

### Simulation flow

`main.rs` loads a JSON case into `Universe` (`src/particles/universe.rs`) and runs an integrator (`src/integrator/whfast.rs` — the standard choice, also `ias15.rs`, `leapfrog.rs`). Each force evaluation calls `Universe::calculate_additional_effects(dangular_momentum_dt, accelerations, ...)`, which sequences per-effect phases: coordinate setup → component precomputation → accelerations → torques (`dangular_momentum_dt`, integrated by a midpoint scheme for spins). Positions are maintained in both inertial and heliocentric frames; each effect keeps its own copy of heliocentric coordinates relative to *its* host particle (`particle.<effect>.coordinates`).

Units: AU, day, solar mass; `constants::K2` is G in these units. `MAX_PARTICLES = 10` (fixed-size arrays throughout; also used to key pair maps as `id * MAX_PARTICLES + depends_on_id`).

### Effects system

Each effect (`tides`, `rotational_flattening`, `general_relativity`, `disk`, `wind`, `evolution` in `src/effects/`) declares per-particle roles via an enum (e.g. `TidesEffect::CentralBody(model) | OrbitingBody(model) | Disabled`). `Universe.hosts` resolves which particle is the host for each effect. `EvolutionType` variants interpolate radius/rg²/k2/1-Q from the `input/` tables each timestep (`src/effects/evolution.rs`).

### Tides (the most intricate effect — `src/effects/tides.rs` + `src/effects/tides/`)

Three tidal models: `ConstantTimeLag` (CTL, Bolmont et al. 2015), `Kaula` (love-number spectra, Revol et al. 2024), `CreepCoplanar` (Gomes et al. 2021).

**Force assembly** (`tides::calculate_tidal_acceleration`) is per-pair and per-body: for each orbiting body, the *planetary-tide* force is dispatched on the orbiting body's model and the *stellar-tide* force on the central body's model, both expressed as force-on-the-planet; the host receives the Newton's-third-law reaction. Any model mixture works, including companions with `TidesEffect::Disabled` (tides coordinates are maintained for **all** particles precisely so the stellar tide survives a tide-disabled companion). Convention trap: `kaula::calculate_tidal_force(a, b, central_body)` swaps roles — with `central_body=true`, argument 1 is the *planet* and the returned vector is the force **on the star** (negated by the caller).

**Per-pair state** lives on the orbiting particle in `tides.parameters.internal`: split orthogonal/radial components `_due_to_stellar_tide` / `_due_to_planetary_tide` (written by CTL and Kaula, consumed by force, torque and `denergy_dt`), and `stellar_tide_secular_force` (Kaula stellar torque, stored per companion — do not collapse it back to a single slot on the star; multi-planet stellar torque depends on it).

**CTL dissipation** σ is composed per `TideComposition` (`Equilibrium | Dynamical | Both`, serde default `Both` = historical behavior): the equilibrium part is `dissipation_factor_scale * dissipation_factor`; the frequency-averaged dynamical part uses the star's evolving `lag_angle` (set from the 1/Q column of BolmontMathis2016 / GalletBolmont2017 / LeconteChabrier2013(true) tables in `evolution.rs`) and only inside the excitation regime `|Ω − n| < Ω`. Active overrides live in `Universe.pair_dependent_scaled_dissipation_factor` (keyed directionally: "dissipation of body `id` as modified by `depends_on_id`").

**Kaula** interpolates complex love numbers linearly over tidal excitation frequencies `w_2mpq = (2−2p+q)n − mΩ` from a 1024-point spectrum (`kaula/love_number.rs`), with a per-timestep cache over (m,p,q) modes and a spin-ratio rescaling for stellar spectra. The q-summation range adapts to eccentricity (`select_eccentricty_order_q`).

Known remaining issues and their history are documented in `docs/known_bugs_ctl_evolution.md`.

### Python side mirrors

`posidonius/effects/tides.py` etc. are thin dict builders mirroring the Rust structs; `posidonius/particles/particle.py` composes them and holds cross-effect consistency warnings (e.g. CTL `tide_composition="Dynamical"` without a 1/Q-providing evolution model). Analysis of history binaries: `posidonius.analysis.history.read(filename) -> (n_particles, records)`; higher-level scripts in `scripts/` (`explore_history.py`, `raw_history.py`, resonance tools).

## Numerical expectations

Rust f64 gives ~16 significant decimals; regression comparisons should use relative tolerances (~1e-12 allows float re-association from refactors while catching physics changes — this is what the golden tests use). Physics-preserving refactors must keep the golden tests green; intended physics changes require regenerating golden values (print helpers are in `golden_tests.rs`) and justifying the change.
