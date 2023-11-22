#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt
from opengate.tests import utility

if __name__ == "__main__":
    """
    This test tets the dose to water functionality relative to NIST values using Silicon as the stopping power ratio and mass stopping power are quite distinct.
    Basic concept: two Water phantoms are calculated. One pure water, one Si with water insert. We compare the ratio of the two phantoms.
    """
    paths = utility.get_default_test_paths(
        __file__, "gate_test041_dose_actor_dose_to_water"
    )

    # create the simulation
    sim = gate.Simulation()

    # main options
    ui = sim.user_info
    ui.g4_verbose = False
    ui.g4_verbose_level = 1
    ui.visu = False
    ui.random_seed = 123456
    ui.number_of_threads = 5
    # units
    # m = gate.g4_units("m")
    # cm = gate.g4_units("cm")
    # mm = gate.g4_units("mm")
    # km = gate.g4_units("km")
    MeV = gate.g4_units.MeV
    Bq = gate.g4_units.Bq
    km = gate.g4_units.km
    cm = gate.g4_units.cm
    mm = gate.g4_units.mm
    gcm3 = gate.g4_units.g_cm3
    kBq = 1000 * Bq

    # add a material database
    # sim.add_material_database(paths.gate_data / "HFMaterials2014.db")

    #  change world size
    world = sim.world
    world.size = [600 * cm, 500 * cm, 500 * cm]
    # world.material = "Vacuum"

    # waterbox
    phantom = sim.add_volume("Box", "phantom")
    phantom.size = [10 * cm, 10 * cm, 10 * cm]
    phantom.translation = [-5 * cm, 0, 0]
    phantom.material = "G4_WATER"
    phantom.color = [0, 0, 1, 1]

    test_material_name = "G4_Si"
    phantom_off = sim.add_volume("Box", "phantom_off")
    phantom_off.mother = phantom.name
    phantom_off.size = [100 * mm, 20 * mm, 20 * mm]
    phantom_off.translation = [0 * mm, 0 * mm, 0 * mm]
    phantom_off.material = test_material_name
    phantom_off.color = [0, 0, 1, 1]

    # water slab
    water_slab_insert = sim.add_volume("Box", "water_slab_insert")
    water_slab_insert.mother = phantom_off.name
    water_slab_insert.size = [2 * mm, 20 * mm, 20 * mm]
    water_slab_insert.translation = [43 * mm, 0, 0]
    water_slab_insert.material = "G4_WATER"
    water_slab_insert.color = [0, 0, 1, 1]
    # si entrance
    Si_entrance_region = sim.add_volume("Box", "Si_entrance_region")
    Si_entrance_region.mother = phantom_off.name
    Si_entrance_region.size = [5 * mm, 20 * mm, 20 * mm]
    Si_entrance_region.translation = [47.5 * mm, 0, 0]
    Si_entrance_region.material = "G4_Si"
    Si_entrance_region.color = [0, 0, 1, 1]

    # physics
    p = sim.get_physics_user_info()
    p.physics_list_name = "QGSP_BIC_EMY"
    sim.physics_manager.global_production_cuts.all = 1000 * km
    # sim.set_cut("world", "all", 1000 * km)

    # default source for tests
    source = sim.add_source("GenericSource", "mysource")
    source.energy.mono = 40 * MeV
    source.particle = "proton"
    source.position.type = "disc"  # pos = Beam, shape = circle + sigma
    # rotate the disc, equiv to : rot1 0 1 0 and rot2 0 0 1
    source.position.rotation = Rotation.from_euler("y", 90, degrees=True).as_matrix()
    # source.position.radius = 8 * mm
    source.position.sigma_x = 2 * mm
    source.position.sigma_y = 2 * mm
    source.position.translation = [0, 0, 0]
    source.direction.type = "momentum"
    source.direction.momentum = [-1, 0, 0]
    source.n = 100

    dose_size = [1000, 1, 1]
    dose_spacing = [0.1, 20.0, 20.0]
    doseActorName_IDD_d = "IDD_d"
    doseActor = sim.add_actor("DoseActor", doseActorName_IDD_d)
    doseActor.output = paths.output / ("test041-" + doseActorName_IDD_d + ".mhd")
    doseActor.mother = phantom_off.name
    doseActor.size = dose_size
    doseActor.spacing = dose_spacing
    doseActor.hit_type = "random"
    doseActor.dose = True

    doseActorName_IDD_d2w = "IDD_d2w"
    doseActorDerived = sim.add_actor("DoseActor", doseActorName_IDD_d2w)
    doseActorDerived.output = paths.output / (
        "test041-" + doseActorName_IDD_d2w + ".mhd"
    )
    doseActorDerived.mother = phantom_off.name
    doseActorDerived.size = doseActor.size
    doseActorDerived.spacing = doseActor.spacing
    doseActorDerived.hit_type = "random"
    doseActorDerived.to_water = True
    doseActorDerived.dose = True

    doseActorName_water_slab_insert_d = "IDD_waterSlab_d"
    doseActorDerived = sim.add_actor("DoseActor", doseActorName_water_slab_insert_d)
    doseActorDerived.output = paths.output / (
        "test041-" + doseActorName_water_slab_insert_d + ".mhd"
    )
    doseActorDerived.mother = water_slab_insert.name
    doseActorDerived.size = doseActor.size
    doseActorDerived.spacing = doseActor.spacing
    doseActorDerived.hit_type = "random"
    doseActorDerived.dose = True

    doseActorName_water_slab_insert_d2w = "IDD_waterSlab_d2w"
    doseActorDerived = sim.add_actor("DoseActor", doseActorName_water_slab_insert_d2w)
    doseActorDerived.output = paths.output / (
        "test041-" + doseActorName_water_slab_insert_d2w + ".mhd"
    )
    doseActorDerived.mother = water_slab_insert.name
    doseActorDerived.size = doseActor.size
    doseActorDerived.spacing = doseActor.spacing
    doseActorDerived.hit_type = "random"
    doseActorDerived.to_water = True
    doseActorDerived.dose = True

    doseActorName_Si_entrance_regiont_d = "IDD_Si_entrance_region_d"
    doseActorDerived = sim.add_actor("DoseActor", doseActorName_Si_entrance_regiont_d)
    doseActorDerived.output = paths.output / (
        "test041-" + doseActorName_Si_entrance_regiont_d + ".mhd"
    )
    doseActorDerived.mother = Si_entrance_region.name
    doseActorDerived.size = doseActor.size
    doseActorDerived.spacing = doseActor.spacing
    doseActorDerived.hit_type = "random"
    doseActorDerived.dose = True

    doseActorName_Si_entrance_regiont_d2w = "IDD_Si_entrance_region_d2w"
    doseActorDerived = sim.add_actor("DoseActor", doseActorName_Si_entrance_regiont_d2w)
    doseActorDerived.output = paths.output / (
        "test041-" + doseActorName_Si_entrance_regiont_d2w + ".mhd"
    )
    doseActorDerived.mother = Si_entrance_region.name
    doseActorDerived.size = doseActor.size
    doseActorDerived.spacing = doseActor.spacing
    doseActorDerived.hit_type = "random"
    doseActorDerived.to_water = True
    doseActorDerived.dose = True

    edepActorName_Si_entrance_regiont_d = "Edep_Si_entrance_region_d"
    edepActor_Si = sim.add_actor("DoseActor", edepActorName_Si_entrance_regiont_d)
    edepActor_Si.output = paths.output / (
        "test041-" + edepActorName_Si_entrance_regiont_d + ".mhd"
    )
    edepActor_Si.mother = Si_entrance_region.name
    edepActor_Si.size = doseActor.size
    edepActor_Si.spacing = doseActor.spacing
    edepActor_Si.hit_type = "random"
    edepActor_Si.dose = False

    edepActorName_Si_entrance_regiont_d2w = "Edep_Si_entrance_region_d2w"
    edepActor_Si_to_water = sim.add_actor(
        "DoseActor", edepActorName_Si_entrance_regiont_d2w
    )
    edepActor_Si_to_water.output = paths.output / (
        "test041-" + edepActorName_Si_entrance_regiont_d2w + ".mhd"
    )
    edepActor_Si_to_water.mother = Si_entrance_region.name
    edepActor_Si_to_water.size = doseActor.size
    edepActor_Si_to_water.spacing = doseActor.spacing
    edepActor_Si_to_water.hit_type = "random"
    edepActor_Si_to_water.to_water = True
    edepActor_Si_to_water.dose = False

    # add stat actor
    s = sim.add_actor("SimulationStatisticsActor", "stats")
    s.track_types_flag = True

    # start simulation
    sim.n = 10
    # output = sim.run()
    output = sim.run(start_new_process=True)

    # print results at the end
    stat = sim.output.get_actor("stats")
    print(stat)

    # ----------------------------------------------------------------------------------------------------------------
    # tests
    print()
    """
    doseFpath_IDD_d = str(
        sim.output.get_actor(doseActorName_IDD_d).user_info.output
    ).replace(".mhd", "-Dose.mhd")
    doseFpath_IDD_d2w = str(
        sim.output.get_actor(doseActorName_IDD_d2w).user_info.output
    ).replace(".mhd", "-Dosetowater.mhd")
    doseFpath_geoWater_d = str(
        sim.output.get_actor(doseActorName_water_slab_insert_d).user_info.output
    ).replace(".mhd", "-Dose.mhd")
    doseFpath_geoWater_d2w = str(
        sim.output.get_actor(doseActorName_water_slab_insert_d2w).user_info.output
    ).replace(".mhd", "-Dosetowater.mhd")

    doseFpath_geoSi_d = str(
        sim.output.get_actor(doseActorName_Si_entrance_regiont_d).user_info.output
    ).replace(".mhd", "-Dose.mhd")
    doseFpath_geoSi_d2w = str(
        sim.output.get_actor(doseActorName_Si_entrance_regiont_d2w).user_info.output
    ).replace(".mhd", "-Dosetowater.mhd")
    """
    doseFpath_IDD_d = sim.output.get_actor(doseActorName_IDD_d).user_info.output
    doseFpath_IDD_d2w = sim.output.get_actor(doseActorName_IDD_d2w).user_info.output
    doseFpath_geoWater_d = sim.output.get_actor(
        doseActorName_water_slab_insert_d
    ).user_info.output
    doseFpath_geoWater_d2w = sim.output.get_actor(
        doseActorName_water_slab_insert_d2w
    ).user_info.output
    doseFpath_geoSi_d = sim.output.get_actor(
        doseActorName_Si_entrance_regiont_d
    ).user_info.output
    doseFpath_geoSi_d2w = sim.output.get_actor(
        doseActorName_Si_entrance_regiont_d2w
    ).user_info.output

    unused = utility.assert_images(
        doseFpath_IDD_d,
        doseFpath_IDD_d2w,
        stat,
        tolerance=100,
        ignore_value=0,
        axis="x",
    )
    mass_ratio_Si_Water = 2.33  # g/cm3
    mass_stopping_power_Water_40MeV = 14.8  # MeV*cm2/g
    mass_stopping_power_Si_40MeV = 11.7  # MeV*cm2/g
    mSPR_40MeV = mass_stopping_power_Water_40MeV / mass_stopping_power_Si_40MeV
    mSPR_40MeV = 1.268771331  # MeV/mm from PSTAR NIST tables, Feb 2023
    mSPR_80MeV = 1.253197674  # from PSTAR NIST tables, Feb 2023
    stopping_power_ratio_Water_Si = mSPR_40MeV / mass_ratio_Si_Water  # 0.542900114

    # utility.warning("Test ratio: dose / dose_to_water in geometry with material: G4_WATER")
    is_ok = utility.assert_images_ratio(
        1.00, doseFpath_geoWater_d, doseFpath_geoWater_d2w, abs_tolerance=0.05
    )

    # utility.warning("Test ratio: dose / dose_to_water in geometry with material: G4_Si")
    is_ok = (
        utility.assert_images_ratio(
            mSPR_40MeV, doseFpath_geoSi_d, doseFpath_geoSi_d2w, abs_tolerance=0.05
        )
        and is_ok
    )

    Fpath_geoSi_edep = sim.output.get_actor(
        edepActorName_Si_entrance_regiont_d
    ).user_info.output
    Fpath_geoSi_edep2water = sim.output.get_actor(
        edepActorName_Si_entrance_regiont_d2w
    ).user_info.output
    is_ok = (
        utility.assert_images_ratio(
            stopping_power_ratio_Water_Si,
            Fpath_geoSi_edep,
            Fpath_geoSi_edep2water,
            abs_tolerance=0.05,
        )
        and is_ok
    )

    utility.test_ok(is_ok)
