import os
import inspect
import shutil
import filecmp
import posidonius
from posidonius.tests import common
from posidonius.tests.test_base import TestBase
import difflib
import json

class Tides_kaula(TestBase):

    def setUp(self):
        TestBase.setUp(self)
        self.current_filename, ignore = os.path.splitext(os.path.basename(__file__)) # Filename without extension
        self.current_dirname = os.path.dirname(inspect.getfile(inspect.currentframe()))

    def tearDown(self):
        TestBase.tearDown(self)

    @staticmethod
    def normalize_json(data):
        """Normalize JSON by converting -0.0 to 0.0 and processing nested structures."""
        if isinstance(data, dict):
            return {k: Tides_kaula.normalize_json(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [Tides_kaula.normalize_json(v) for v in data]
        elif isinstance(data, float) and data == -0.0:
            return 0.0
        else:
            return data

    def test_enabled_star_tides(self):

        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_effects = posidonius.ConsiderEffects({
            "tides": True,
            "rotational_flattening": False,
            "general_relativity": False,
            "disk": False,
            "wind": False,
            "evolution": False,
        })
        # general_relativity_implementation = "Kidder1995" # Assumes one central massive body
        #general_relativity_implementation = "Anderson1975" # Assumes one central massive body
        #general_relativity_implementation = "Newhall1983" # Considers all bodies
        general_relativity_implementation = "Disabled" # Assumes one central massive body
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

        star_mass = 1. # Solar masses
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        #star_evolution = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        # star_evolution = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        # star_evolution = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        star_evolution = posidonius.NonEvolving()
        # star = common.brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, general_relativity_implementation, star_evolution)
        star = common.solar_like_for_kaula_1(star_mass, star_position, star_velocity, general_relativity_implementation, star_evolution)
        universe.add_particle(star)

        common.basic_configuration_kaula_1(universe)

        ############################################################################
        # whfast_alternative_coordinates="Jacobi"
        whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        # universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        # Load and normalize both JSON files
        with open(json_filename, 'r') as file1, open(expected_json_filename, 'r') as file2:
            json1 = self.normalize_json(json.load(file1))
            json2 = self.normalize_json(json.load(file2))

        # Compare normalized JSON and generate a diff
        are_files_equal = json1 == json2

        if not are_files_equal:
            diff = difflib.unified_diff(
                json.dumps(json1, indent=4).splitlines(),
                json.dumps(json2, indent=4).splitlines(),
                fromfile=json_filename,
                tofile=expected_json_filename,
                lineterm=''
            )
            print("Differences between JSON files:")
            print('\n'.join(diff))

        # Assert that the JSON files are equal
        self.assertTrue(are_files_equal, "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))

    def test_enabled_planet_tides(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_effects = posidonius.ConsiderEffects({
            "tides": True,
            "rotational_flattening": False,
            "general_relativity": False,
            "disk": False,
            "wind": False,
            "evolution": False,
        })
        # general_relativity_implementation = "Kidder1995" # Assumes one central massive body
        #general_relativity_implementation = "Anderson1975" # Assumes one central massive body
        #general_relativity_implementation = "Newhall1983" # Considers all bodies
        general_relativity_implementation = "Disabled" # Assumes one central massive body
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

        star_mass = 0.0898 # Solar masses
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        #star_evolution = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        # star_evolution = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        # star_evolution = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        star_evolution = posidonius.NonEvolving()
        # star = common.brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, general_relativity_implementation, star_evolution)
        star = common.red_dwarf_like_for_kaula(star_mass, star_position, star_velocity, general_relativity_implementation, star_evolution)
        universe.add_particle(star)

        common.basic_configuration_kaula_2(universe)

        ############################################################################
        # whfast_alternative_coordinates="Jacobi"
        whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        # universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        # Load and normalize both JSON files
        with open(json_filename, 'r') as file1, open(expected_json_filename, 'r') as file2:
            json1 = self.normalize_json(json.load(file1))
            json2 = self.normalize_json(json.load(file2))

        # Compare normalized JSON and generate a diff
        are_files_equal = json1 == json2

        if not are_files_equal:
            diff = difflib.unified_diff(
                json.dumps(json1, indent=4).splitlines(),
                json.dumps(json2, indent=4).splitlines(),
                fromfile=json_filename,
                tofile=expected_json_filename,
                lineterm=''
            )
            print("Differences between JSON files:")
            print('\n'.join(diff))

        # Assert that the JSON files are equal
        self.assertTrue(are_files_equal, "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))

        # Clean up temporary files
        shutil.rmtree(os.path.dirname(json_filename))


    def test_enabled_both_tides(self):
        current_function_name = inspect.stack()[0][3][5:] # Remove first 5 characters "test_"
        expected_json_filename, json_filename = common.setup(self.current_dirname, self.current_filename, current_function_name)

        initial_time, time_step, time_limit, historic_snapshot_period, recovery_snapshot_period = common.simulation_properties()
        consider_effects = posidonius.ConsiderEffects({
            "tides": True,
            "rotational_flattening": False,
            "general_relativity": False,
            "disk": False,
            "wind": False,
            "evolution": False,
        })
        # general_relativity_implementation = "Kidder1995" # Assumes one central massive body
        #general_relativity_implementation = "Anderson1975" # Assumes one central massive body
        #general_relativity_implementation = "Newhall1983" # Considers all bodies
        general_relativity_implementation = "Disabled" # Assumes one central massive body
        universe = posidonius.Universe(initial_time, time_limit, time_step, recovery_snapshot_period, historic_snapshot_period, consider_effects)

        star_mass = 1. # Solar masses
        star_position = posidonius.Axes(0., 0., 0.)
        star_velocity = posidonius.Axes(0., 0., 0.)

        #star_evolution = posidonius.Baraffe1998(star_mass) # M-Dwarf (mass = 0.10) or SolarLike ConstantDissipation (mass = 1.0)
        #star_evolution = posidonius.Baraffe2015(star_mass) # mass = 0.01 .. 1.40
        # star_evolution = posidonius.Leconte2011(star_mass) # BrownDwarf (mass = 0.01 .. 0.08)
        # star_evolution = posidonius.BolmontMathis2016(star_mass) # SolarLike Evolving dissipation (mass = 0.40 .. 1.40)
        #star_evolution = posidonius.GalletBolmont2017(star_mass) # SolarLike Evolving dissipation (mass = 0.30 .. 1.40)
        #star_evolution = posidonius.LeconteChabrier2013(False) # Jupiter without dissipation of dynamical tides
        #star_evolution = posidonius.LeconteChabrier2013(True) # Jupiter with dissipation of dynamical tides
        star_evolution = posidonius.NonEvolving()
        # star = common.brown_dwarf(star_mass, star_dissipation_factor_scale, star_position, star_velocity, general_relativity_implementation, star_evolution)
        star = common.solar_like_for_kaula_2(star_mass, star_position, star_velocity, general_relativity_implementation, star_evolution)
        universe.add_particle(star)

        common.basic_configuration_kaula_2(universe)

        ############################################################################
        # whfast_alternative_coordinates="Jacobi"
        whfast_alternative_coordinates="DemocraticHeliocentric"
        #whfast_alternative_coordinates="WHDS"
        #universe.write(json_filename, integrator="LeapFrog")
        # universe.write(json_filename, integrator="IAS15")
        universe.write(json_filename, integrator="WHFast", whfast_alternative_coordinates=whfast_alternative_coordinates)

        if not os.path.exists(expected_json_filename):
            shutil.copyfile(json_filename, expected_json_filename)

        # Load and normalize both JSON files
        with open(json_filename, 'r') as file1, open(expected_json_filename, 'r') as file2:
            json1 = self.normalize_json(json.load(file1))
            json2 = self.normalize_json(json.load(file2))

        # Compare normalized JSON and generate a diff
        are_files_equal = json1 == json2

        if not are_files_equal:
            diff = difflib.unified_diff(
                json.dumps(json1, indent=4).splitlines(),
                json.dumps(json2, indent=4).splitlines(),
                fromfile=json_filename,
                tofile=expected_json_filename,
                lineterm=''
            )
            print("Differences between JSON files:")
            print('\n'.join(diff))

        # Assert that the JSON files are equal
        self.assertTrue(are_files_equal, "Generated JSON case is not equal to the expected one: {} {}".format(json_filename, expected_json_filename))

        # Clean up temporary files
        shutil.rmtree(os.path.dirname(json_filename))

