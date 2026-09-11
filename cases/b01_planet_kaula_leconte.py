"""Planetary Kaula tide only, mirrored in Spiroid (see ~/Documents/spi_pos_comparison/run_spiroid_b01.py).
Non-evolving 1 Msun star (same kg mass as Spiroid's SOLAR_MASS), Venus-mass planet at 0.03 AU, e = 0.05, coplanar,
planet spin = 1.03 x mean motion, Leconte+2015/Steinberger k2 spectrum read from the Spiroid JSON so both codes use
identical data. Usage: python3 cases/b01_planet_kaula_leconte.py out.json [--notides] [--days 300]"""
import posidonius, numpy as np, argparse, json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename')
    parser.add_argument('--notides', action='store_true')
    parser.add_argument('--days', type=float, default=300.)
    parser.add_argument('--spectrum', default='../spiroid/examples/data/planet/tides/kaula/leconte2015_steinberger.json')
    args = parser.parse_args()
    C = posidonius.constants
    tides_on = not args.notides

    initial_time = 5e6 * 365.25; time_step = 0.01; time_limit = args.days
    universe = posidonius.Universe(initial_time, time_limit, time_step, 365.25 * 1e5, 0.5,
                                   posidonius.ConsiderEffects({"tides": tides_on, "rotational_flattening": False, "general_relativity": False,
                                                               "disk": False, "wind": False, "evolution": False}))
    # --- star: same kg mass as Spiroid (astro-const SOLAR_MASS = 1.32712440099e20 / 6.67428e-11)
    star_mass_kg = 1.32712440099e20 / 6.67428e-11
    star_mass = star_mass_kg / C.M_SUN
    star_spin_z = C.TWO_PI / 1.91                       # days^-1 (irrelevant for the planetary tide)
    star = posidonius.Particle(star_mass, 1.0 * C.R_SUN, 0.2645751311, posidonius.Axes(0., 0., 0.), posidonius.Axes(0., 0., 0.), posidonius.Axes(0., 0., star_spin_z))
    star.set_tides(posidonius.effects.tides.CentralBody(posidonius.effects.tides.DisabledModel()))
    star.set_rotational_flattening(posidonius.effects.rotational_flattening.Disabled())
    star.set_general_relativity(posidonius.effects.general_relativity.Disabled())
    star.set_wind(posidonius.effects.wind.Disabled())
    star.set_disk(posidonius.effects.disk.Disabled())
    star.set_evolution(posidonius.NonEvolving())
    universe.add_particle(star)
    # --- planet (Venus mass/radius/rg2 as in spiroid/examples/planet_kaula_solid_tides_1d_interpolation.conf)
    planet_mass_kg = 4.8685e24; planet_radius_m = 6052000.0; rg2 = 0.33070368308499226
    planet_mass = planet_mass_kg / C.M_SUN; planet_radius = planet_radius_m / C.AU
    a = 0.03; e = 0.05; i = 0.; p = 0.; n = 0.; l = 0.
    q = a * (1.0 - e)
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star.get_position()] if hasattr(star, 'get_position') else [posidonius.Axes(0., 0., 0.)], velocities=[posidonius.Axes(0., 0., 0.)])
    mean_motion = np.sqrt(C.G_SI * (star_mass_kg + planet_mass_kg) / (a * C.AU)**3)   # rad/s
    planet_spin_rad_s = 1.03 * mean_motion
    planet_spin = posidonius.Axes(0., 0., planet_spin_rad_s * C.DAY)                  # days^-1
    # --- k2 spectrum: Spiroid JSON (x_vals rad/s, data [Re, Im]); pad to 1024 beyond the last frequency with the edge value
    sp = json.load(open(args.spectrum))["Interpolate1D"]
    w = np.array(sp["x_vals"], float); k2 = np.array(sp["data"], float); Re = k2[:, 0]; Im = k2[:, 1]
    npad = 1024 - len(w); assert npad >= 0
    w = np.concatenate([w, w[-1] + (w[-1] - w[-2]) * np.arange(1, npad + 1)]); Re = np.concatenate([Re, [Re[-1]] * npad]); Im = np.concatenate([Im, [Im[-1]] * npad])
    assert np.all(np.diff(w) > 0)
    kaula = posidonius.effects.tides.Kaula({"love_numbers": {"spectrum_excitation_frequency": w.tolist(), "spectrum_real_part": Re.tolist(), "spectrum_imaginary_part": Im.tolist()}})
    planet = posidonius.Particle(planet_mass, planet_radius, np.sqrt(rg2), planet_position, planet_velocity, planet_spin)
    planet.set_tides(posidonius.effects.tides.OrbitingBody(kaula if tides_on else posidonius.effects.tides.DisabledModel()))
    planet.set_rotational_flattening(posidonius.effects.rotational_flattening.Disabled())
    planet.set_general_relativity(posidonius.effects.general_relativity.Disabled())
    planet.set_wind(posidonius.effects.wind.Disabled())
    planet.set_disk(posidonius.effects.disk.Disabled())
    planet.set_evolution(posidonius.NonEvolving())
    universe.add_particle(planet)
    print("star mass %.10e kg, planet mass %.6e kg, radius %.6e m, a %.6e m, e %.3f, n %.10e rad/s, planet spin %.10e rad/s" % (star_mass_kg, planet_mass_kg, planet_radius_m, a * C.AU, e, mean_motion, planet_spin_rad_s))
    universe.write(args.output_filename, integrator="WHFast", whfast_alternative_coordinates="DemocraticHeliocentric")
