import numpy as np
from collections import namedtuple
from scipy.spatial.transform import Rotation as R
import opengate as gate
from ...definitions import __world_name__
from ...element import copy_user_info

Solid = namedtuple(
    "Solid", "solid_type properties"
)  # TODO: in principle we should use the solid constructor
Volume = namedtuple("Volume", "material solid")


def gdml_to_gate_volume_name(gdml_name):
    return gdml_name.capitalize()  # + "Volume"


def get_list_from_xml_element(element, unit_name="unit"):
    unit = gate.g4_units[element.get(unit_name)]
    vec = []
    for key, value in element.items():
        if key in ["name", unit_name]:
            continue
        else:
            vec.append(float(value) * unit)

    return vec


def get_solids_gdml(root, positions=None):
    # print(f'--- In {element.tag} ---')
    solids = {}
    for solids_el in root.findall("solids"):
        for child in solids_el:
            solid_type = gdml_to_gate_volume_name(child.tag)
            solid_name = child.get("name")
            if solid_type == "Tessellated":
                vectors = []
                for triangle in child.findall("triangular"):
                    vectors.append(get_tessellated_solid_vector(positions, triangle))

                properties = vectors

            else:
                properties = get_list_from_xml_element(child, unit_name="lunit")

            solids[solid_name] = Solid(solid_type, properties)

    return solids


def get_volumes_gdml(root, solids):
    volumes = {}
    for structure in root.findall("structure"):
        for volume in structure.findall("volume"):
            volume_name = volume.get("name")
            material = volume.find("materialref").get("ref")
            solidref = volume.find("solidref").get("ref")
            volumes[volume_name] = Volume(material, solids[solidref])
    return volumes


def add_volumes_gdml(vol_manager, root):
    positions = get_positions_gdml(root)
    solids = get_solids_gdml(root, positions=positions)
    volumes = get_volumes_gdml(root, solids)
    assemblies = get_assemblies(root)
    world_name = get_world_name(root)

    # create a volume for each physical volume
    for phv in assemblies[world_name]:
        add_volume(
            phv, volumes, vol_manager, root, assemblies, mother_name=__world_name__
        )


def add_volume(phv, volumes, vol_manager, root, assemblies, mother_name):
    name = phv.get("name")
    volumeref = phv.find("volumeref").get("ref")
    if volumeref in volumes:
        if name in vol_manager.volume_names:
            name += "1"  # repetition of same solid
        volume = volumes[volumeref]
        volume_type = volume.solid.solid_type
        # print(name, mother_name)
        # add new volume
        vol = vol_manager.add_volume(volume_type, name)
        if volume_type == "Tessellated":
            vol.vectors = volume.solid.properties
        else:
            vol.size = volume.solid.properties
        # add material
        vol.material = volume.material
        # add rotation and translation if present
        vol.translation = get_translation_gdml(phv)
        vol.rotation = get_rotation_gdml(phv)
        # hierarchy
        vol.mother = mother_name
    elif volumeref in assemblies:
        vol = vol_manager.add_volume("Box", name)
        # print(name, get_translation_gdml(phv))
        vol.translation = get_translation_gdml(phv)
        vol.rotation = get_rotation_gdml(phv)
        for phv in assemblies[volumeref]:
            add_volume(phv, volumes, vol_manager, root, assemblies, mother_name=name)


def get_world_name(root):
    for setup in root.findall("setup"):
        world_name = setup.find("world").get("ref")
    return world_name


def get_assemblies(root):
    assemblies = {}
    for structure in root.findall("structure"):
        for assembly in structure.findall("assembly"):
            assembly_name = assembly.get("name")
            phys_vols = assembly.findall("physvol")
            assemblies[assembly_name] = phys_vols
    return assemblies


def get_rotation_gdml(phv):
    # rotation is in  Tait-Bryan angles (pitch, yaw, roll)
    r = phv.find("rotation")
    rotation = R.identity().as_matrix()

    # OSS: element will test false also if it doesn't have sub elements! (if p: syntax is deprecated)
    if r is not None:
        rot_vec = get_list_from_xml_element(r)
        deg = False if r.get("unit") == "rad" else True
        rotation = R.from_euler("yzx", rot_vec, degrees=deg).as_matrix()
    return rotation


def get_translation_gdml(phv):
    p = phv.find("position")
    translation = [0.0, 0.0, 0.0]
    # OSS: element will test false also if it doesn't have sub elements! (if p: syntax is deprecated)
    if p is not None:
        # print(phv.get('name'))
        translation = get_list_from_xml_element(p)

    return translation


def get_tessellated_solid_vector(positions, triangle):
    v1 = positions[triangle.get("vertex1")]
    v2 = positions[triangle.get("vertex2")]
    v3 = positions[triangle.get("vertex3")]
    vec = [v1, v2, v3]

    return np.array(vec)


def get_positions_gdml(root):
    positions = {}
    for p in root.iter("position"):
        positions[p.get("name")] = np.array(get_list_from_xml_element(p))

    return positions
