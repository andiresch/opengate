#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt
import numpy as np


def run_sim(N_events: int, N_threads: int, N_voxels: int, paths):
    # create the simulation
    sim = gate.Simulation()

    # main options
    ui = sim.user_info
    ui.g4_verbose = False
    ui.g4_verbose_level = 1
    ui.visu = False
    # ui.random_seed = 123456789
    ui.number_of_threads = N_threads
    Ntotal = N_events
    N_per_trhead = int(np.round(Ntotal / ui.number_of_threads))
    # units
    m = gate.g4_units("m")
    cm = gate.g4_units("cm")
    mm = gate.g4_units("mm")
    km = gate.g4_units("km")
    MeV = gate.g4_units("MeV")
    Bq = gate.g4_units("Bq")
    kBq = 1000 * Bq

    # add a material database
    # sim.add_material_database(paths.gate_data / "HFMaterials2014.db")

    #  change world size
    world = sim.world
    world.size = [600 * cm, 500 * cm, 500 * cm]
    # world.material = "Vacuum"

    # waterbox
    phantom = sim.add_volume("Box", "phantom")
    phantom.size = [50 * mm, 100 * mm, 100 * mm]
    phantom.translation = [-25 * mm, 0, 0]
    phantom.material = "G4_WATER"
    phantom.color = [0, 0, 1, 1]

    # physics
    p = sim.get_physics_user_info()
    p.physics_list_name = "QGSP_BIC_EMY"
    sim.physics_manager.global_production_cuts.all = 1000 * km

    # default source for tests
    source = sim.add_source("GenericSource", "mysource")
    source.energy.mono = 80 * MeV
    source.particle = "proton"
    source.position.type = "disc"
    source.position.rotation = Rotation.from_euler("y", 90, degrees=True).as_matrix()
    source.position.sigma_x = 2 * mm
    source.position.sigma_y = 2 * mm
    source.position.translation = [0, 0, 0]
    source.direction.type = "momentum"
    source.direction.momentum = [-1, 0, 0]
    source.n = N_per_trhead

    # define actors:
    dose_size = [N_voxels, 1, 1]
    # dose_spacing = [5 * mm, 100.0 * mm, 100.0 * mm]
    dose_spacing = np.divide(phantom.size, dose_size)
    print(dose_spacing)

    # create default/reference image
    doseActorName_IDD_singleImage = "IDD_singleImage"
    doseActor = sim.add_actor("DoseActor", doseActorName_IDD_singleImage)
    doseActor.output = paths.output / (
        "test041-" + doseActorName_IDD_singleImage + ".mhd"
    )
    doseActor.mother = phantom.name
    doseActor.size = dose_size
    doseActor.spacing = dose_spacing
    doseActor.hit_type = "random"
    doseActor.dose = False
    doseActor.use_more_RAM = False
    doseActor.ste_of_mean = False
    doseActor.uncertainty = False
    doseActor.square = False

    # create test image: img for each thread
    doseActorName_IDD_NthreadImages = "IDD_NthreadImages"
    doseActor = sim.add_actor("DoseActor", doseActorName_IDD_NthreadImages)
    doseActor.output = paths.output / (
        "test041-" + doseActorName_IDD_NthreadImages + ".mhd"
    )
    doseActor.mother = phantom.name
    doseActor.size = dose_size
    doseActor.spacing = dose_spacing
    doseActor.hit_type = "random"
    doseActor.dose = False
    doseActor.use_more_RAM = True
    doseActor.ste_of_mean = False
    doseActor.uncertainty = False
    doseActor.square = False

    # add stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True

    # start simulation
    sim.n = int(N_per_trhead)

    output = sim.run(start_new_process=True)

    # print results at the end
    stat = sim.output.get_actor("stats")
    print(stat)

    # ----------------------------------------------------------------------------------------------------------------
    # tests
    doseFpath_IDD_singleImage = sim.output.get_actor(
        doseActorName_IDD_singleImage
    ).user_info.output
    doseFpath_IDD_NthreadImages = sim.output.get_actor(
        doseActorName_IDD_NthreadImages
    ).user_info.output
    return doseFpath_IDD_singleImage, doseFpath_IDD_NthreadImages, stat


def run_test(doseFpath_IDD_singleImage, doseFpath_IDD_NthreadImages, stat):
    unused = gate.assert_images(
        doseFpath_IDD_singleImage,
        doseFpath_IDD_NthreadImages,
        stat,
        tolerance=100,
        ignore_value=0,
        axis="x",
    )
    expected_ratio = 1.00
    gate.warning("Test ratio: edep in single image vs. Nthread image")
    is_ok = gate.assert_images_ratio(
        expected_ratio,
        doseFpath_IDD_singleImage,
        doseFpath_IDD_NthreadImages,
        abs_tolerance=0.03,
    )
    return is_ok


if __name__ == "__main__":
    paths = gate.get_default_test_paths(
        __file__, "gate_test041_dose_actor_dose_to_water"
    )
    N_threadsV = [8, 2]
    N_voxelsV = [50, 2]
    N_events = 800
    summary_str = []
    is_ok = True
    out_str = ""
    for N_thread in N_threadsV:
        for N_voxel in N_voxelsV:
            fPathRef_singleImage, fPathTest_NThreadImgs, stat = run_sim(
                N_events=N_events, N_threads=N_thread, N_voxels=N_voxel, paths=paths
            )
            is_ok_i = run_test(fPathRef_singleImage, fPathTest_NThreadImgs, stat)
            out_str += f"Test {is_ok_i}, for {N_thread = } and {N_voxel = }" + "\n"

            is_ok = is_ok_i and is_ok

    print(out_str)

    gate.test_ok(is_ok)
