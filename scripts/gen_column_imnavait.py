#!/usr/bin/env python3
"""
Generate a 1D column mesh for Imnavait Creek with vertical stratification.

This script creates a vertical mesh with layered organic and mineral soils,
starting from the surface elevation (892 m). The mesh includes:
    - Acrotelm (6 cm)
    - Catotelm (27 cm)
    - Shallow and deep mineral soil (with expanding dz)
    - Reaches to a total depth of 40 m

The mesh is generated as a single-cell surface (1D column) and can be used
in ATS simulations. The vertical discretization is exported to .npy for 
reference or use in analysis.

Environment variables required:
    - AMANZI_TPLS_DIR
    - ATS_SRC_DIR

"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

# Load Exodus (SEACAS) module
try:
    import exodus
except ImportError:
    sys.path.append(os.path.join(os.environ['AMANZI_TPLS_DIR'], 'SEACAS', 'lib'))
    import exodus

# Load ATS meshing tools
try:
    import meshing_ats
except ImportError:
    sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'meshing', 'meshing_ats'))
    import meshing_ats

# ----------------------------
# SURFACE MESH: Single Column
# ----------------------------

# Flat surface at elevation 892 m
x = np.array([0.0, 1.0], dtype='d')       # Horizontal extent (1 m wide)
z = np.array([892.00, 892.00], dtype='d') # Constant surface elevation

# Construct 2D transect (actually 1D column, single cell in y-dir)
m2 = meshing_ats.Mesh2D.from_Transect(x, z)

# ----------------------------
# LAYER SETUP
# ----------------------------

layer_types = []
layer_data = []
layer_ncells = []
layer_mat_ids = []
z_curr = 0.0

print("Layer | Cell ID | Mat ID | dz [m]")

# --- Acrotelm (6 cm, 1 cm resolution)
layer_types.append("constant")
layer_data.append(0.06)
layer_ncells.append(6)
layer_mat_ids.append(1001)
z_curr -= 0.06

# Print layer info
count = 0
for i, thick in enumerate(layer_data):
    for j in range(layer_ncells[i]):
        print(f"{i:>5} | {count:>7} | {layer_mat_ids[i]:>6} | {thick/layer_ncells[i]:.4f}")
        count += 1

# --- Catotelm (27 cm, 1 cm resolution)
layer_types.append("constant")
layer_data.append(0.27)
layer_ncells.append(27)
layer_mat_ids.append(1002)
z_curr -= 0.27

# --- Mineral Layer Part 1 (60 cm, 1 cm resolution)
layer_types.append("constant")
layer_data.append(0.6)
layer_ncells.append(60)
layer_mat_ids.append(1003)
z_curr -= 0.6

# --- Mineral Layer Part 2: 50 cells with increasing dz
dz = 0.02
for _ in range(50):
    dz *= 1.05
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1003)
    z_curr -= dz

# --- Mineral Layer Part 3: Continue increasing until dz reaches ~1.8 m
while dz < 1.8:
    dz *= 1.05
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1003)
    z_curr -= dz

# Print new layers
for i in range(3, len(layer_data)):  # Skip first 3 already printed
    for j in range(layer_ncells[i]):
        print(f"{i:>5} | {count:>7} | {layer_mat_ids[i]:>6} | {layer_data[i]/layer_ncells[i]:.4f}")
        count += 1

# --- Mineral Layer Part 4: Constant 2.0 m cells until 40 m total depth
remaining_depth = 40 + z_curr  # z_curr is negative
layer_types.append("constant")
layer_data.append(remaining_depth)
layer_ncells.append(int(round(remaining_depth / 2.0)))
layer_mat_ids.append(1003)

# ----------------------------
# Save vertical grid resolution
# ----------------------------

cell_widths = []
for i, thick in enumerate(layer_data):
    dz = thick / layer_ncells[i]
    cell_widths.extend([dz] * layer_ncells[i])

cell_widths = np.array(cell_widths).reshape(-1, 1)
np.save("cell_widths_v2.npy", cell_widths)

# ----------------------------
# EXTRUDE 3D MESH AND SAVE
# ----------------------------

m3 = meshing_ats.Mesh3D.extruded_Mesh2D(
    m2,
    layer_types,
    layer_data,
    layer_ncells,
    layer_mat_ids
)

# Uncomment to save the mesh file

m3.write_exodus("../data/mesh_column_imnavait.exo")
print("Exodus mesh written to: ../data/mesh_column_imnavait.exo")
print(f"\nTotal number of vertical cells: {cell_widths.shape[0]}")
print("Vertical grid saved to: cell_widths_v2.npy")