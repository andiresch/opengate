import opengate as gate
import itk
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

## ------ INITIALIZE SIMULATION ENVIRONMENT ---------- ##
paths = gate.get_default_test_paths(__file__, "gate_test044_pbs")
output_path = paths.output / "output_test051_rtp"
ref_path = paths.output_ref / "test051_ref"


# create the simulation
sim = gate.Simulation()

# main options
ui = sim.user_info
ui.g4_verbose = False
ui.g4_verbose_level = 1
ui.visu = False
ui.random_seed = 12365478910
ui.random_engine = "MersenneTwister"

# units
km = gate.g4_units("km")
cm = gate.g4_units("cm")
mm = gate.g4_units("mm")

# add a material database
sim.add_material_database(paths.gate_data / "HFMaterials2014.db")

#  change world size
world = sim.world
world.size = [600 * cm, 500 * cm, 500 * cm]

## ---------- DEFINE BEAMLINE MODEL -------------##
beamline = gate.BeamlineModel()
beamline.Name = None
beamline.RadiationTypes = "ion 6 12"
# Nozzle entrance to Isocenter distance
beamline.NozzleToIsoDist = 1300.00  # 1648 * mm#1300 * mm
# SMX to Isocenter distance
beamline.SMXToIso = 6700.00
# SMY to Isocenter distance
beamline.SMYToIso = 7420.00
# polinomial coefficients
beamline.energyMeanCoeffs = [-6.71618e-9, 1.02304e-5, -4.9270e-3, 1.28461e1, -66.136]
beamline.energySpreadCoeffs = [-1.66295e-9, 1.31502e-6, -2.59769e-4, -2.60088e-3, 7.436]
beamline.sigmaXCoeffs = [
    -1.07268e-13,
    1.61558e-10,
    -9.92211e-8,
    3.19029e-5,
    -5.67757e-3,
    5.29884e-1,
    -17.749,
]
beamline.thetaXCoeffs = [
    -1.13854e-17,
    1.52020e-14,
    -7.49359e-12,
    1.57991e-9,
    -8.98373e-8,
    -1.30862e-5,
    1.638e-3,
]
beamline.epsilonXCoeffs = [
    -2.54669e-16,
    3.71028e-13,
    -2.14188e-10,
    6.21900e-8,
    -9.46711e-6,
    7.09187e-4,
    -19.511e-3,
]
beamline.sigmaYCoeffs = [
    -5.80689e-14,
    9.10249e-11,
    -5.75230e-8,
    1.85977e-5,
    -3.20430e-3,
    2.74490e-1,
    -7.133,
]
beamline.thetaYCoeffs = [
    8.10201e-18,
    -1.75709e-14,
    1.44445e-11,
    -5.82592e-9,
    1.22471e-6,
    -1.28547e-4,
    6.066e-3,
]
beamline.epsilonYCoeffs = [
    -5.74235e-16,
    9.12245e-13,
    -5.88501e-10,
    1.96763e-7,
    -3.58265e-5,
    3.35307e-3,
    -122.935e-3,
]

# NOTE: HBL means that the beam is coming from -x (90 degree rot around y)
nSim = 20000  # particles to simulate per beam
# rt_plan = ref_path / "RP1.2.752.243.1.1.20220406175810679.4500.52008_tagman.dcm"
# beamset = gate.beamset_info(rt_plan)
# G = float(beamset.beam_angles[0])


## ----  VBL Nozzle  ---
# nozzle box
box = sim.add_volume("Box", "box")
box.size = [500 * mm, 500 * mm, 1000 * mm]
box.translation = [0.0 * mm, 0.0, 1148.0]
box.material = "Vacuum"
box.color = [0, 0, 1, 1]

# nozzle WET
nozzle = sim.add_volume("Box", "nozzle")
nozzle.mother = box.name
nozzle.size = [500 * mm, 500 * mm, 2 * mm]
nozzle.material = "G4_WATER"

# Rashi
rashi = sim.add_volume("Box", "rashi")
rashi.mother = box.name
rashi.size = [500 * mm, 500 * mm, 5 * mm]
rashi.translation = [0.0, 0.0, -200 * mm]
rashi.material = "G4_LUCITE"
rashi.color = [1, 0, 1, 1]

## ----  HBL Nozzle  ---
box_rot = sim.add_volume("Box", "box_rot")
gate.copy_user_info(box, box_rot)
box_rot.rotation = Rotation.from_euler("y", 90, degrees=True).as_matrix()
box_rot.translation = [1148.0, 1000.0, 0.0]

nozzle_rot = sim.add_volume("Box", "nozzle_rot")
gate.copy_user_info(nozzle, nozzle_rot)
nozzle_rot.mother = box_rot.name

rashi_rot = sim.add_volume("Box", "rashi_rot")
gate.copy_user_info(rashi, rashi_rot)
rashi_rot.mother = box_rot.name

# -----------------------------------

# target 1
phantom = sim.add_volume("Box", "phantom")
phantom.size = [324 * mm, 324 * mm, 324 * mm]
phantom.translation = [0 * mm, 0.0, 0.0]
phantom.material = "G4_WATER"
phantom.color = [0, 0, 1, 1]

# target 2
phantom_rot = sim.add_volume("Box", "phantom_rot")
gate.copy_user_info(phantom, phantom_rot)
phantom_rot.translation = [0.0, 1000.0, 0.0]
phantom_rot.rotation = Rotation.from_euler("y", 90, degrees=True).as_matrix()

# add dose actor
dose = sim.add_actor("DoseActor", "doseInXYZ")
dose.output = paths.output / "testTPSxyz.mhd"
dose.mother = phantom.name
dose.size = [162, 162, 648]
dose.spacing = [2.0, 2.0, 2.0]
dose.hit_type = "random"
dose.gray = True

dose_rot = sim.add_actor("DoseActor", "doseInXYZ_rot")
gate.copy_user_info(dose, dose_rot)
dose_rot.mother = phantom_rot.name
dose_rot.output = paths.output / "testTPSxyz_rot.mhd"

# physics
p = sim.get_physics_user_info()
p.physics_list_name = "FTFP_INCLXX_EMZ"
sim.set_cut("world", "all", 1000 * km)

# add TPSources
spots, ntot, energies, G = gate.spots_info_from_txt(
    ref_path / "TreatmentPlan4Gate-1D_HBL_120.txt", "ion 6 12"
)
tps = gate.TreatmentPlanSource(nSim, sim, beamline)
# tps.beamset = beamset
tps.spots = spots
tps.name = "VBL"
tps.rotation = Rotation.from_euler("y", 0, degrees=True)
tps.initialize_tpsource()

tps_rot = gate.TreatmentPlanSource(nSim, sim, beamline)
tps_rot.spots = spots
tps_rot.name = "HBL"
tps_rot.rotation = Rotation.from_euler("y", G, degrees=True)
tps_rot.translation = [0.0, 1000.0, 0.0]
tps_rot.initialize_tpsource()

# add stat actor
s = sim.add_actor("SimulationStatisticsActor", "Stats")
s.track_types_flag = True
# start simulation
output = sim.start()

# print results at the end
stat = output.get_actor("Stats")
print(stat)

# create output dir, if it doesn't exist
if not os.path.isdir(output_path):
    os.mkdir(output_path)

## ------ TESTS -------##

# ABSOLUTE DOSE
ok = gate.assert_images(
    dose.output,
    dose_rot.output,
    stat,
    tolerance=50,
    ignore_value=0,
)

# read output and ref
img_mhd_out = itk.imread(dose_rot.output)
img_mhd_ref = itk.imread(dose.output)
data = itk.GetArrayViewFromImage(img_mhd_out)
data_ref = itk.GetArrayViewFromImage(img_mhd_ref)
spacing = img_mhd_out.GetSpacing()

# Range 80
ok = (
    gate.compareRange(data, data_ref, data.shape, data_ref.shape, spacing, spacing)
    and ok
)

# 1D plots
fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 10))
gate.plot_img_axis(ax, img_mhd_out, "z profile", axis="z")
gate.plot_img_axis(ax, img_mhd_out, "x profile", axis="x")
gate.plot_img_axis(ax, img_mhd_out, "y profile", axis="y")

# fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(25, 10))
gate.plot_img_axis(ax, img_mhd_ref, "z ref", axis="z")
gate.plot_img_axis(ax, img_mhd_ref, "x ref", axis="x")
gate.plot_img_axis(ax, img_mhd_ref, "y ref", axis="y")
# fig.savefig(output_path / "dose_profiles_water.png")
# plt.show()

gate.test_ok(ok)