#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
from opengate.tests import utility
from scipy.spatial.transform import Rotation
import numpy as np
import uproot


# The test generates two differents generic sources, momentum and focused, which is attached to a plan.
# The plan is randomly rotated and we verify that the generated particles have a direction which is in accordance
# with the applied rotations and the transmissions.


def test068(tab, nb_run, nb_part, nb_source):
    nb_event = len(tab)
    nb_event_theo = nb_run * nb_part * nb_source
    err = np.sqrt(nb_event_theo)
    if nb_event_theo - 3 * err < nb_event < nb_event_theo + 3 * err:
        return True
    else:
        return False


if __name__ == "__main__":
    paths = utility.get_default_test_paths(__file__)
    output_path = paths.output

    # create the simulation
    sim = gate.Simulation()

    # main options
    ui = sim.user_info
    ui.g4_verbose = False
    # ui.visu = True
    ui.visu_type = "vrml"
    ui.check_volumes_overlap = False
    # ui.running_verbose_level = gate.logger.EVENT
    ui.number_of_threads = 1
    ui.random_seed = "auto"

    # units
    m = gate.g4_units.m
    km = gate.g4_units.km
    mm = gate.g4_units.mm
    cm = gate.g4_units.cm
    nm = gate.g4_units.nm
    Bq = gate.g4_units.Bq
    MeV = gate.g4_units.MeV
    keV = gate.g4_units.keV
    sec = gate.g4_units.s
    gcm3 = gate.g4_units["g/cm3"]

    #  adapt world size
    world = sim.world
    world.size = [1 * m, 1 * m, 1 * m]
    world.material = "G4_Galactic"

    # nb of particles
    n = 100

    # Plan to attach the source
    nb_source = 2
    plan = sim.add_volume("Box", "source_plan")
    plan.mother = world.name
    plan.material = "G4_Galactic"
    plan.size = [10 * cm, 10 * cm, 1 * nm]
    plan.translation = np.array([0 * cm, 0 * cm, 0 * cm])

    source1 = sim.add_source("GenericSource", "photon_source")
    source1.particle = "gamma"
    source1.position.type = "box"
    source1.mother = plan.name
    source1.position.size = [10 * cm, 10 * cm, 1 * nm]
    source1.direction.type = "momentum"
    source1.force_rotation = True
    # source1.direction.focus_point = [0*cm, 0*cm, -5 *cm]
    source1.direction.momentum = [0, 0, -1]
    source1.energy.type = "mono"
    source1.energy.mono = 1 * MeV
    source1.activity = n * Bq / ui.number_of_threads

    source2 = sim.add_source("GenericSource", "photon_source_2")
    source2.particle = "gamma"
    source2.position.type = "box"
    source2.mother = plan.name
    source2.position.size = [10 * cm, 10 * cm, 1 * nm]
    source2.direction.type = "focused"
    source2.force_rotation = True
    source2.direction.focus_point = [0 * cm, 0 * cm, -5 * cm]
    # source1.direction.momentum = [0, 0, -1]
    source2.energy.type = "mono"
    source2.energy.mono = 1 * MeV
    source2.activity = n * Bq / ui.number_of_threads

    # Phase Space

    phsp_plan = sim.add_volume("Box", "phsp_plan")
    phsp_plan.mother = world.name
    phsp_plan.material = "G4_Galactic"
    phsp_plan.size = [10 * cm, 10 * cm, 1 * nm]
    phsp_plan.translation = np.array([0 * cm, 0 * cm, -5 * cm])

    phsp = sim.add_actor("PhaseSpaceActor", "PhaseSpace")
    phsp.mother = phsp_plan.name
    phsp.attributes = ["EventID"]
    phsp.output = output_path / "test068.root"

    # MOTION ACTOR for the source and actors
    motion_source = sim.add_actor("MotionVolumeActor", "Move_source")
    motion_source.mother = plan.name
    motion_source.rotations = []
    motion_source.translations = []

    motion_phsp = sim.add_actor("MotionVolumeActor", "Move_phsp")
    motion_phsp.mother = phsp_plan.name
    motion_phsp.rotations = []
    motion_phsp.translations = []
    sim.run_timing_intervals = []

    for i in range(10):
        x_translation = 20 * cm * np.random.random()
        y_translation = 20 * cm * np.random.random()
        z_translation = 20 * cm * np.random.random()

        x_rot = 360 * np.random.random()
        y_rot = 360 * np.random.random()
        z_rot = 360 * np.random.random()

        motion_source.translations.append([x_translation, y_translation, z_translation])
        rot = Rotation.from_euler(
            "xyz", [x_rot, y_rot, z_rot], degrees=True
        ).as_matrix()
        vec_source_phsp = phsp_plan.translation - plan.translation
        # vec_source_phsp = vec_source_phsp.reshape((len(vec_source_phsp),1))
        # print(vec_source_phsp)
        translation_phsp = np.dot(rot, vec_source_phsp) + np.array(
            [x_translation, y_translation, z_translation]
        )
        motion_source.rotations.append(rot)
        motion_phsp.rotations.append(rot)
        motion_phsp.translations.append(translation_phsp)
        sim.run_timing_intervals.append([i * sec, (i + 1) * sec])
    # stat actor
    s = sim.add_actor("SimulationStatisticsActor", "Stats")
    s.track_types_flag = True

    # Physic list and cuts
    sim.physics_manager.physics_list_name = "G4EmStandardPhysics_option3"
    sim.physics_manager.enable_decay = False
    sim.physics_manager.global_production_cuts.gamma = 1 * mm
    sim.physics_manager.global_production_cuts.electron = 1 * mm
    sim.physics_manager.global_production_cuts.positron = 1 * mm

    # go !
    sim.run()
    output = sim.output

    # Validation test
    f_phsp = uproot.open(phsp.output)
    arr = f_phsp["PhaseSpace"].arrays()
    print("Number of detected events :", len(arr))
    nb_event = len(arr)
    nb_part = n
    nb_run = len(sim.run_timing_intervals)
    err = np.sqrt(nb_part * nb_run * nb_source)
    print(
        "Number of expected events : ["
        + str(int(-3 * err + nb_run * nb_part * nb_source))
        + ","
        + str(int(+3 * err + nb_run * nb_part * nb_source))
        + "]"
    )
    is_ok = test068(arr, len(sim.run_timing_intervals), nb_part, 2)
    utility.test_ok(is_ok)
