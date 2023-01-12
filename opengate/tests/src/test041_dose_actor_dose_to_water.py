#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import opengate as gate
from scipy.spatial.transform import Rotation
import matplotlib.pyplot as plt

paths = gate.get_default_test_paths(__file__, "gate_test042_gauss_gps")

# create the simulation
sim = gate.Simulation()

# main options
ui = sim.user_info
ui.g4_verbose = False
ui.g4_verbose_level = 1
ui.visu = False
ui.random_seed = 123456

# units
m = gate.g4_units("m")
cm = gate.g4_units("cm")
mm = gate.g4_units("mm")
km = gate.g4_units("km")
MeV = gate.g4_units("MeV")
Bq = gate.g4_units("Bq")
kBq = 1000 * Bq

# add a material database
sim.add_material_database(paths.gate_data / "HFMaterials2014.db")

#  change world size
world = sim.world
world.size = [600 * cm, 500 * cm, 500 * cm]
world.material = "Vacuum"

# waterbox
phantom = sim.add_volume("Box", "phantom")
phantom.size = [10 * cm, 20 * cm, 20 * cm]
phantom.translation = [-5 * cm, 0, 0]
phantom.material = "G4_WATER"
phantom.color = [0, 0, 1, 1]

# daughter
phantom_y = sim.add_volume("Box", "phantom_y")
phantom_y.mother = phantom.name
phantom_y.size = [100 * mm, 20 * mm, 20 * mm]
# phantom_y.translation = [49 * mm, 0, 0]
phantom_y.material = "G4_WATER"
phantom_y.color = [0, 0, 1, 1]

test_material_name = "G4_Si"
phantom_off = sim.add_volume("Box", "phantom_off")
phantom_off.mother = phantom.name
phantom_off.size = [100 * mm, 20 * mm, 20 * mm]
phantom_off.translation = [0 * mm, 50 * mm, 0 * mm]
phantom_off.material = test_material_name
phantom_off.color = [0, 0, 1, 1]


# daughter
slab_insert = sim.add_volume("Box", "slab_insert")
slab_insert.mother = phantom_off.name
slab_insert.size = [5 * mm, 20 * mm, 20 * mm]
slab_insert.translation = [40 * mm, 0, 0]
slab_insert.material = "G4_WATER"
slab_insert.color = [0, 0, 1, 1]


# physics
p = sim.get_physics_user_info()
p.physics_list_name = "QGSP_BIC_EMY"
sim.set_cut("world", "all", 1000 * km)
# FIXME need SetMaxStepSizeInRegion ActivateStepLimiter

# default source for tests
source = sim.add_source("Generic", "mysource")
source.energy.mono = 80 * MeV
source.particle = "proton"
source.position.type = "disc"  # pos = Beam, shape = circle + sigma
# rotate the disc, equiv to : rot1 0 1 0 and rot2 0 0 1
source.position.rotation = Rotation.from_euler("y", 90, degrees=True).as_matrix()
# source.position.radius = 8 * mm
source.position.sigma_x = 8 * mm
source.position.sigma_y = 4 * mm
source.position.translation = [0, 50, 0]
source.direction.type = "momentum"
source.direction.momentum = [-1, 0, 0]
source.activity = 100 * kBq

# add dose actor
dose = sim.add_actor("DoseActor", "doseInXZ")
dose.output = paths.output / "test041-lateral_xz.mhd"
dose.mother = phantom_y.name
dose.size = [250, 1, 250]
dose.spacing = [0.4, 100, 0.4]
dose.hit_type = "random"
#

dose = sim.add_actor("DoseActor", "doseInXY")
dose.output = paths.output / "test041-lateral_xy.mhd"
dose.mother = phantom_y.name
dose.size = [250, 250, 1]
dose.spacing = [0.4, 0.4, 100]
dose.hit_type = "random"

doseActorName3 = "doseInYZ"
dose = sim.add_actor("DoseActor", doseActorName3)
dose.output = paths.output / ("test041-" + doseActorName3 + ".mhd")
dose.mother = phantom_y.name
# dose.size = [1, 250, 250]
# dose.spacing = [100, 0.4, 0.4]
dose.size = [100, 1, 1]
dose.spacing = [1.0, 20.0, 20.0]
dose.hit_type = "random"


doseActorName4 = "doseActor4"
doseFour = sim.add_actor("DoseActor", doseActorName4)
doseFour.output = paths.output / ("test041-" + doseActorName4 + ".mhd")
doseFour.mother = phantom_off.name
# dose.size = [1, 250, 250]
# dose.spacing = [100, 0.4, 0.4]
doseFour.size = [100, 1, 1]
doseFour.spacing = [1.0, 20.0, 20.0]
doseFour.hit_type = "random"
doseFour.gray = True


doseActorName5 = "doseActor5"
doseFive = sim.add_actor("DoseActor", doseActorName5)
doseFive.output = paths.output / ("test041-" + doseActorName5 + ".mhd")
doseFive.mother = phantom_off.name
# dose.size = [1, 250, 250]
# dose.spacing = [100, 0.4, 0.4]
doseFive.size = [100, 1, 1]
doseFive.spacing = [1.0, 20.0, 20.0]
doseFive.hit_type = "random"
doseFive.gray = True
doseFive.dose_to_water = True

# add stat actor
s = sim.add_actor("SimulationStatisticsActor", "stats")
s.track_types_flag = True

# start simulation
sim.n = 10
output = sim.start()

# print results at the end
stat = output.get_actor("stats")
print(stat)

# dose = output.get_actor("doseInXZ")
# print(dose)
# plt.plot(dose)
# plt.show()

# ----------------------------------------------------------------------------------------------------------------
# tests
print()
gate.warning("Tests stats file")
stats_ref = gate.read_stat_file(paths.gate_output / "stats.txt")
is_ok = gate.assert_stats(stat, stats_ref, 0.14)

spr_m = 1.18
# is_ok = gate.test_weights(
#    spr_m, output_path / mhd_1, output_path / mhd_2, thresh=0.2
# )
"""
print('hier:::::::::::', output.get_actor("doseInXZ").user_info.output)
gate.warning("Difference for EDEP XZ")
is_ok = (
    gate.assert_images(
        output.get_actor("doseInXZ").user_info.output,
        #paths.gate_output / "lateral_xz_Protons_40MeV_sourceShapeGaussian-Edep.mhd",
        output.get_actor("doseInXZ").user_info.output,
        stat,
        tolerance=10,
        ignore_value=0,
    )
    and is_ok
)

print()
gate.warning("Difference for EDEP XY")
is_ok = (
    gate.assert_images(
        output.get_actor("doseInXY").user_info.output,
        output.get_actor("doseInXY").user_info.output,
        stat,
        tolerance=10,
        ignore_value=0,
        axis="y",
    )
    and is_ok
)

print()
"""
gate.warning("Difference for EDEP YZ")
fName4_dose = str(output.get_actor(doseActorName4).user_info.output).replace(
    ".mhd", "_dose.mhd"
)
fName5_dosewater = str(output.get_actor(doseActorName5).user_info.output).replace(
    ".mhd", "_doseToWater.mhd"
)
is_ok = (
    gate.assert_images(
        fName4_dose,
        fName5_dosewater,
        stat,
        tolerance=30,
        ignore_value=0,
        axis="x",
    )
    and is_ok
)
mSPR_40MeV = 1.268771331
mSPR_80MeV = 1.253197674

is_ok = gate.test_weights(1.26, fName4_dose, fName5_dosewater, thresh=0.2)


gate.test_ok(is_ok)
