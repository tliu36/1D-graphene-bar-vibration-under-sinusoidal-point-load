# -*- coding: utf-8 -*-
from odbAccess import *
import numpy as np

# === Input Parameters ===
odb_path = '20250814_bar_load_1kHz_f_1_625n5.odb'
instance_name = 'BAR-1'
step_name = 'Step-1'
output_file = '1_20250814_bar_load_1kHz_f_1_625n5_u2_0_16000.txt'

start_frame = 0
end_frame = 10000  # inclusive
z_tol = 1e-4  

# === Open the ODB ===
odb = openOdb(path=odb_path, readOnly=True)
step = odb.steps[step_name]
instance = odb.rootAssembly.instances[instance_name]

# === Get Z-coordinates of all nodes ===
node_z_coords = {}
for node in instance.nodes:
    node_z_coords[node.label] = node.coordinates[2]  # Z is 3rd coordinate

# === Open output file for writing ===
with open(output_file, 'w') as f:

    for frame_id in range(start_frame, end_frame + 1):
        frame = step.frames[frame_id]
        displacement_field = frame.fieldOutputs['U']
        disp_subset = displacement_field.getSubset(region=instance, position=NODAL)

        f.write("Frame {}\n".format(frame_id))

        seen_nodes = set()
        data = []

        for val in disp_subset.values:
            node_id = val.nodeLabel
            if node_id in seen_nodes:
                continue
            seen_nodes.add(node_id)

            u2 = val.data[1]  # Y-direction
            z = node_z_coords.get(node_id, None)
            if z is not None:
                data.append((node_id, z, u2))

        # === Sort by node ID to match original logic ===
        data.sort(key=lambda x: x[0])

        # === Average every 2 lines and write only the average ===
        for i in range(0, len(data)-1, 2):
            id_new = i//2 + 1
            _, z1, u2_1 = data[i]
            _, z2, u2_2 = data[i + 1]
            z_avg = 0.5 * (z1 + z2)
            u2_avg = 0.5 * (u2_1 + u2_2)
            f.write("{:d} {:.6f} {:.6e}\n".format(id_new, z_avg, u2_avg))

# === Close the ODB ===
odb.close()