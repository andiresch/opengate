import numpy as np


def gdml_to_gate_volume_name(gdml_name):
    return gdml_name.capitalize() + "Volume"


def initialize_solids_gdml(sim_engine, element, positions=None):
    # print(f'--- In {element.tag} ---')
    for child in element:
        # print(f'  --- In {child.tag} ---')
        # print(child.tag, child.attrib)
        if element.tag == "solids":
            volume_type = gdml_to_gate_volume_name(child.tag)
            volume_name = child.get("name")
            vol = sim_engine.add_volume(volume_type, volume_name)
            if volume_type == "TessellatedVolume":
                vectors = []
                for triangle in child.findall("triangular"):
                    vectors.append(get_tessellated_solid_vector(positions, triangle))

                vol.vectors = vectors

            else:
                size = [
                    float(child.get("x")),
                    float(child.get("y")),
                    float(child.get("z")),
                ]
                print(volume_type, volume_name, size)
                vol.size = size

        initialize_solids_gdml(sim_engine, child, positions=positions)


def get_volumes_gdml(sim_engine, root):
    for child in root:
        for volume in child.findall("volume"):
            material = volume.find("materialref").get("ref")
            solidref = volume.find("solidref").get("ref")
            print(solidref, material)

            if solidref in sim_engine.volume_names:
                sim_engine.volumes[solidref].material = material


def get_tessellated_solid_vector(positions, triangle):
    v1 = np.array(positions[triangle.get("vertex1")])
    v2 = np.array(positions[triangle.get("vertex2")])
    v3 = np.array(positions[triangle.get("vertex3")])
    vec = [v1, v2, v3]

    return np.array(vec)


def get_positions_gdml(root):
    positions = {}
    for p in root.iter("position"):
        positions[p.get("name")] = np.array(
            [float(p.get("x")), float(p.get("y")), float(p.get("z"))]
        )

    return positions
