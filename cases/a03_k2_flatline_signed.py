import posidonius
import numpy as np
import argparse
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('output_filename', action='store', help='Filename where the initial snapshot will be stored (e.g., universe_integrator.json)')
    parser.add_argument('--spectrum', action='store', default='../spiroid/examples/data/star/tides/spectrum_1msun_10myr_flat.json',
                        help='Stellar k2 spectrum JSON in the Spiroid format (default: the flat spectrum shared with Spiroid)')
    parser.add_argument('--evolution', choices=['galletbolmont', 'starevol'], default='galletbolmont',
                        help='stellar evolution track: GalletBolmont2017 (default) or the Starevol table shared with Spiroid (input/Starevol/M_10.dat)')
    parser.add_argument('--days', type=float, default=3., help='time limit in days (default 3)')
    parser.add_argument('--snapshot', type=float, default=None, help='historic snapshot period in days (default: the time step)')
    parser.add_argument('--ecc', type=float, default=0., help='initial eccentricity (default 0)')
    parser.add_argument('--notides', action='store_true', help='disable the stellar tide (control run)')
    parser.add_argument('--second_planet_au', type=float, default=0., help='add a distant Earth-mass planet (tides disabled) at this semi-major axis in AU (0 = none); used to test the multi-planet stellar torque')

    args = parser.parse_args()
    filename = args.output_filename

    initial_time = 5e6*365.25 # time [days] where simulation starts
    time_step = 0.01 # days
    time_limit = args.days  # days
    historic_snapshot_period = args.snapshot if args.snapshot is not None else time_step # days
    recovery_snapshot_period = 365.25 * 100000. # days
    consider_effects = posidonius.ConsiderEffects({
        "tides": not args.notides,
        "rotational_flattening": False,
        "general_relativity": False,
        "disk": False,
        "wind": True,
        "evolution": True,
    })
    universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

    # Solar masses. With the Starevol track the star mass is set to Spiroid's SOLAR_MASS (1.32712440099e20 / G) in kg
    # so that both codes use the same numbers (Spiroid overrides the conf mass with table_mass x SOLAR_MASS).
    star_mass = (1.32712440099e20 / 6.67428e-11) / posidonius.constants.M_SUN if args.evolution == 'starevol' else 1.
    star_radius_factor = 1.#0.117
    star_radius = star_radius_factor * posidonius.constants.R_SUN
    star_radius_of_gyration = 0.2645751311; # Sun, so alpha = 0.07
    star_position = posidonius.Axes(0., 0., 0.)
    star_velocity = posidonius.Axes(0., 0., 0.)

    # Initialization of stellar spin
    star_rotation_period = 1.91*24 # hours (original 3.3)
    star_angular_frequency = posidonius.constants.TWO_PI/(star_rotation_period/24.) # days^-1
    star_spin = posidonius.Axes(0., 0., star_angular_frequency)

    # Stellar k2 spectrum: read the SAME file Spiroid uses (spiroid/examples/data/star/tides/),
    # so that both codes see an identical love number spectrum.
    #   - The JSON stores dimensionless tidal frequencies x = omega / spin_spec, with spin_spec in rad/s.
    #   - Posidonius interpolates in SI (rad/s, see kaula.rs: spin and orbital_frequency are divided by DAY),
    #     so the axis is converted back to rad/s by multiplying with spin_spec.
    #   - spin_spec is also the "stellar_spectrum_spin_rate" used for the (spin / spin_spec)^2 rescaling
    #     of Im(k2), the same rescaling Spiroid applies (love_number.rs, StarConvectiveEnvelope).
    #   - data = [Re(k2), Im(k2)]. Both codes negate Re(k2) internally and use Im(k2) at the SIGNED frequency
    #     for the star (no parity flip), so the values are passed through unchanged.
    with open(args.spectrum) as f:
        spectrum_json = json.load(f)
    spectrum_spin_rate = float(spectrum_json["spin_spec"])               # rad/s
    spectrum_table = spectrum_json["spectrum"]["Interpolate1D"]
    x_spec   = np.array(spectrum_table["x_vals"], dtype=float)           # omega / spin_spec
    k2_spec  = np.array(spectrum_table["data"], dtype=float)             # columns: Re, Im
    ReK2_spec = k2_spec[:, 0]
    ImK2_spec = k2_spec[:, 1]

    # Posidonius expects exactly 1024 points. Resample each sign of the frequency separately so that
    # the sign change of Im(k2) at x = 0 stays as sharp as in the source table (no interpolation across 0).
    expected_size = 1024
    half = expected_size // 2
    neg = x_spec < 0.
    pos = x_spec > 0.
    x_neg = np.linspace(x_spec[neg].min(), x_spec[neg].max(), half)
    x_pos = np.linspace(x_spec[pos].min(), x_spec[pos].max(), half)
    w_lm_star = np.concatenate([x_neg, x_pos]) * spectrum_spin_rate     # rad/s, strictly increasing
    ReK2_star = np.concatenate([np.interp(x_neg, x_spec[neg], ReK2_spec[neg]),
                                np.interp(x_pos, x_spec[pos], ReK2_spec[pos])])
    ImK2_star = np.concatenate([np.interp(x_neg, x_spec[neg], ImK2_spec[neg]),
                                np.interp(x_pos, x_spec[pos], ImK2_spec[pos])])
    size_star = np.size(w_lm_star)
    assert size_star == expected_size
    assert np.all(np.diff(w_lm_star) > 0.)
    print("Spectrum: %s" % args.spectrum)
    print("  spin_spec = %.10e rad/s (P = %.6f d)" % (spectrum_spin_rate, posidonius.constants.TWO_PI / spectrum_spin_rate / posidonius.constants.DAY))
    print("  omega range = [%.4e, %.4e] rad/s, Re(k2) in [%.6g, %.6g], Im(k2) in [%.4e, %.4e]"
          % (w_lm_star.min(), w_lm_star.max(), ReK2_star.min(), ReK2_star.max(), ImK2_star.min(), ImK2_star.max()))

    star_kaula_tidal_parameters_love_numbers = {
        "love_numbers": {
            "spectrum_excitation_frequency": w_lm_star.tolist(),
            "spectrum_real_part": ReK2_star.tolist(),
            "spectrum_imaginary_part": ImK2_star.tolist(),
            "stellar_spectrum_spin_rate": spectrum_spin_rate
        },
    }

    star_tides_model = posidonius.effects.tides.Kaula(star_kaula_tidal_parameters_love_numbers)
    # star_tides_model = posidonius.effects.tides.ConstantTimeLag(star_tides_parameters)
    star_tides = posidonius.effects.tides.CentralBody(star_tides_model if not args.notides else posidonius.effects.tides.DisabledModel())
    # star_tides = posidonius.effects.tides.OrbitingBody(star_tides_model)
    # star_tides = posidonius.effects.tides.Disabled()
    #
    # star_rotational_flattening_parameters = {"love_number": star_tides_parameters["love_number"]}
    # star_rotational_flattening_parameters =  {"uniform_viscosity_coefficient": 1e14}
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(star_rotational_flattening_parameters)
    # star_rotational_flattening_model = posidonius.effects.rotational_flattening.CreepCoplanar(star_rotational_flattening_parameters)
    # star_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(star_rotational_flattening_model)
    #star_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(star_rotational_flattening_model)
    star_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    #
    # star_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #star_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #star_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #star_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    star_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    star_wind = posidonius.effects.wind.Interaction({
        # Solar wind parametrisation (Bouvier 1997)
        "k_factor": 5.15e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
        "rotation_saturation": 9. * posidonius.constants.TWO_PI / 25.0, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    })
#     star_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
    #'inner_edge_distance': 0.01,  # AU
    #'outer_edge_distance': 100.0, # AU
    #'lifetime': 1.0e5 * 365.25e0, # days
    #'alpha': 1.0e-2,
    #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
    #'mean_molecular_weight': 2.4,
    #}
    #star_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #star_disk = posidonius.effects.disk.OrbitingBody()
    star_disk = posidonius.effects.disk.Disabled()
    #
    star_evolution = posidonius.Starevol(star_mass) if args.evolution == 'starevol' else posidonius.GalletBolmont2017(star_mass) # mass = 0.30 .. 1.40
    # star_evolution = posidonius.BolmontMathis2016(star_mass) # mass = 0.40 .. 1.40
    #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
    #star_evolution = posidonius.Leconte2011(star_mass) # mass = 0.01 .. 0.08
    #star_evolution = posidonius.Baraffe1998(star_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
    # star_evolution = posidonius.NonEvolving()
    #
    star = posidonius.Particle(star_mass, star_radius, star_radius_of_gyration, star_position, star_velocity, star_spin)
    star.set_tides(star_tides)
    star.set_rotational_flattening(star_rotational_flattening)
    star.set_general_relativity(star_general_relativity)
    star.set_wind(star_wind)
    star.set_disk(star_disk)
    star.set_evolution(star_evolution)
    universe.add_particle(star)

    ############################################################################

    Mearth = 5.9736E+24 #Earth Mass (kg)
    Rearth = 6378.1E+3    #Earth radius (m)
    Mvenus = 4.8685E+24
    Rvenus = 6052E+3
    Mjup = 1.898E+27 # Jupiter mass (kg)
    Rjup = 6.991E+07 # Jupiter radius (m)

    # Radius computed from equation of Zeng 2016
    planet_radius_of_gyration = 5.75e-01 # Earth type planet
    planet_mass = 5. * posidonius.constants.M_EARTH # Solar masses (3.0e-6 solar masses = 1 earth mass)
    planet_radius = 1.5687291659 * posidonius.constants.R_EARTH # Assume same average density with Earth's

    # first planet
    #////////// Specify initial position and velocity for a stable orbit
    #////// Keplerian orbital elements, in the `asteroidal' format of Mercury code
    a = 0.03                # semi-major axis (in AU)
    e = args.ecc;                               # eccentricity
    i = 0. * posidonius.constants.DEG2RAD;  # inclination (degrees)
    p = 0. * posidonius.constants.DEG2RAD;  # argrument of pericentre (degrees)
    n = 0. * posidonius.constants.DEG2RAD;  # longitude of the ascending node (degrees)
    l = 0. * posidonius.constants.DEG2RAD;  # mean anomaly (degrees)
    p = (p + n);                            # Convert to longitude of perihelion !!
    q = a * (1.0 - e);                      # perihelion distance
    planet_position, planet_velocity = posidonius.calculate_cartesian_coordinates(planet_mass, q, e, i, p, n, l, masses=[star_mass], positions=[star_position], velocities=[star_velocity])

    #////// Initialization of planetary spin
    planet_obliquity = 0. * posidonius.constants.DEG2RAD # 0.2 rad
    orbital_period = 2.*np.pi*(np.sqrt( (a*posidonius.constants.AU)**3. / ( posidonius.constants.G_SI*( star_mass*posidonius.constants.M_SUN + planet_mass*posidonius.constants.M_SUN ) ) ) )/86400.
    planet_rotation_period = orbital_period / (1.03) #24. # hours orbital_period/(1.0 + 1e-5) #24. # hours
    planet_angular_frequency = posidonius.constants.TWO_PI/(planet_rotation_period) # days^-1
    planet_keplerian_orbital_elements = posidonius.calculate_keplerian_orbital_elements(planet_mass, planet_position, planet_velocity, masses=[star_mass], positions=[star_position], velocities=[star_velocity])
    planet_inclination = planet_keplerian_orbital_elements[3]
    (planet_sma, planet_peri_distance, planet_ecc, planet_inclination, planet_long_peri, planet_long_asc_node, planet_mean_anomaly, planet_true_ano) = planet_keplerian_orbital_elements
    planet_spin = posidonius.calculate_spin(planet_angular_frequency, planet_inclination, planet_obliquity, planet_long_asc_node)

    # k2pdelta = 2.465278e-3 # Terrestrial planets (no gas)
    # planet_tides_parameters = {
    #     "dissipation_factor_scale": 1.0,
    #     "dissipation_factor": 2. * posidonius.constants.K2 * k2pdelta/(3. * np.power(planet_radius, 5)),
    #     "love_number": 0.299,
    # }

    # planet_tides_parameters = {
    #     "dissipation_factor_scale": 0.0,
    #     "dissipation_factor": 0.0,
    #     "love_number": 0.0,
    # }

    # --- Choose the tidal model to use :
    # planet_tides_model = posidonius.effects.tides.Kaula(planet_kaula_tidal_parameters_love_numbers)
    # planet_tides_model = posidonius.effects.tides.ConstantTimeLag(planet_tides_parameters)

    # --- Choose the type of particle (central body or orbiting body) :
    #planet_tides = posidonius.effects.tides.CentralBody(planet_tides_model)
    # planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
    planet_tides_model = posidonius.effects.tides.DisabledModel()
    planet_tides = posidonius.effects.tides.OrbitingBody(planet_tides_model)
    # planet_tides = posidonius.effects.tides.Disabled()
    # ---
    #planet_rotational_flattening_parameters = {"love_number": 0.9532}
    #planet_rotational_flattening_model = posidonius.effects.rotational_flattening.OblateSpheroid(planet_rotational_flattening_parameters)
    #planet_rotational_flattening = posidonius.effects.rotational_flattening.CentralBody(planet_rotational_flattening_model)
    # planet_rotational_flattening = posidonius.effects.rotational_flattening.OrbitingBody(planet_rotational_flattening_model)
    planet_rotational_flattening = posidonius.effects.rotational_flattening.Disabled()
    # ---
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Kidder1995")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Anderson1975")
    #planet_general_relativity = posidonius.effects.general_relativity.CentralBody("Newhall1983")
    #planet_general_relativity = posidonius.effects.general_relativity.OrbitingBody()
    planet_general_relativity = posidonius.effects.general_relativity.Disabled()
    #
    #planet_wind = posidonius.effects.wind.Interaction({
    ## Solar wind parametrisation (Bouvier 1997)
    #"k_factor": 4.0e-18, # K_wind = 1.6d47 cgs, which is in Msun.AU2.day
    #"rotation_saturation": 1.7592918860102842, # 14. * TWO_PI/25.0, in units of the spin of the Sun today
    #})
    planet_wind = posidonius.effects.wind.Disabled()
    #
    #disk_surface_density_normalization_gcm = 1000. # g.cm^-2
    #disk_surface_density_normalization_SI = disk_surface_density_normalization_gcm * 1.0e-3 * 1.0e4 # kg.m^-2
    #disk_properties = {
    #'inner_edge_distance': 0.01,  # AU
    #'outer_edge_distance': 100.0, # AU
    #'lifetime': 1.0e5 * 365.25e0, # days
    #'alpha': 1.0e-2,
    #'surface_density_normalization': disk_surface_density_normalization_SI * (1.0/posidonius.constants.M_SUN) * posidonius.constants.AU**2, # Msun.AU^-2
    #'mean_molecular_weight': 2.4,
    #}
    #planet_disk = posidonius.effects.disk.CentralBody(disk_properties)
    #planet_disk = posidonius.effects.disk.OrbitingBody()
    planet_disk = posidonius.effects.disk.Disabled()
    #
    #planet_evolution = posidonius.GalletBolmont2017(planet_mass) # mass = 0.30 .. 1.40
    #planet_evolution = posidonius.BolmontMathis2016(planet_mass) # mass = 0.40 .. 1.40
    #planet_evolution = posidonius.Baraffe2015(planet_mass) # mass = 0.01 .. 1.40
    #planet_evolution = posidonius.Leconte2011(planet_mass) # mass = 0.01 .. 0.08
    #planet_evolution = posidonius.Baraffe1998(planet_mass) # Sun (mass = 1.0) or M-Dwarf (mass = 0.1)
    #planet_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
    #planet_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
    planet_evolution = posidonius.NonEvolving()
    #
    planet = posidonius.Particle(planet_mass, planet_radius, planet_radius_of_gyration, planet_position, planet_velocity, planet_spin)
    planet.set_tides(planet_tides)
    planet.set_rotational_flattening(planet_rotational_flattening)
    planet.set_general_relativity(planet_general_relativity)
    planet.set_wind(planet_wind)
    planet.set_disk(planet_disk)
    planet.set_evolution(planet_evolution)
    universe.add_particle(planet)

    if args.second_planet_au > 0.:
        # Distant Earth-mass planet without tides: its stellar tide is negligible, so the stellar torque must
        # stay that of the first planet (README_debug_kaula_2026.md, F9).
        p2_mass = 1. * posidonius.constants.M_EARTH
        p2_radius = 1. * posidonius.constants.R_EARTH
        p2_position, p2_velocity = posidonius.calculate_cartesian_coordinates(p2_mass, args.second_planet_au, 0., 0., 0., 0., 0., masses=[star_mass], positions=[star_position], velocities=[star_velocity])
        p2_n = posidonius.constants.TWO_PI / (2. * np.pi * np.sqrt((args.second_planet_au * posidonius.constants.AU)**3 / (posidonius.constants.G_SI * star_mass * posidonius.constants.M_SUN)) / 86400.)
        p2 = posidonius.Particle(p2_mass, p2_radius, 0.5, p2_position, p2_velocity, posidonius.Axes(0., 0., p2_n))
        p2.set_tides(posidonius.effects.tides.OrbitingBody(posidonius.effects.tides.DisabledModel()))
        p2.set_rotational_flattening(posidonius.effects.rotational_flattening.Disabled())
        p2.set_general_relativity(posidonius.effects.general_relativity.Disabled())
        p2.set_wind(posidonius.effects.wind.Disabled())
        p2.set_disk(posidonius.effects.disk.Disabled())
        p2.set_evolution(posidonius.NonEvolving())
        universe.add_particle(p2)

    ############################################################################

    whfast_alternative_coordinates="DemocraticHeliocentric"
    #whfast_alternative_coordinates="WHDS"
    #whfast_alternative_coordinates="Jacobi"
    universe.write(filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)
#     universe.write(filename, integrator="IAS15")
    # universe.write(filename, integrator="LeapFrog")
