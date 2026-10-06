#!/usr/bin/env python3
"""
Generate a 3D extruded mesh for a hillslope transect at Imnavait Creek.

This script constructs a 2D transect profile from surveyed x-z coordinates
and extrudes it vertically with stratified soil layers including:
    - Acrotelm
    - Catotelm
    - Organic/mineral layers with constant and expanding vertical resolution

The final mesh is written in Exodus format, ready for use in ATS simulations.

Requirements:
    - meshing_ats
    - SEACAS (Exodus Python bindings)
    - Environment variables: AMANZI_TPLS_DIR, ATS_SRC_DIR
"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt

# -- Load Exodus Python bindings --
try:
    import exodus
except ImportError:
    sys.path.append(os.path.join(os.environ['AMANZI_TPLS_DIR'], 'SEACAS', 'lib'))
    import exodus

# -- Load ATS meshing tools --
try:
    import meshing_ats
except ImportError:
    sys.path.append(os.path.join(os.environ['ATS_SRC_DIR'], 'tools', 'meshing', 'meshing_ats'))
    import meshing_ats

# ----------------------
# 2D TRANSECT GEOMETRY
# ----------------------

# Measured topographic profile (x in meters, z in meters above sea level)
x = np.array([0, 9.33, 22.31, 30.30, 35.87, 41.84, 47.36, 53.44, 59.45,
              65.51, 71.86, 77.70, 84.20, 90.42, 95.91, 98.68, 102.63, 106.68, 109.00])
z = np.array([892.00, 891.71, 891.41, 891.24, 891.17, 891.02, 890.93, 890.81, 890.58,
              890.41, 890.33, 890.17, 889.93, 889.77, 889.58, 889.53, 889.41, 889.25, 889.00])

# Interpolate to regular spacing along x
n = 109  # desired horizontal resolution
x_interp = np.linspace(x[0], x[-1], n)
z_interp = np.interp(x_interp, x, z)
x, z = x_interp, z_interp
dx = x[n-5] - x[n-6]

# Generate 2D mesh
m2 = meshing_ats.Mesh2D.from_Transect(x, z)

# ----------------------
# EXTRUSION LAYERS
# ----------------------

layer_types = []
layer_data = []
layer_ncells = []
layer_mat_ids = []

z_curr = 0.0  # current depth (used only for tracking)

# Acrotelm (0.06 m, 1 cm cells)
layer_types.append("constant")
layer_data.append(0.06)
layer_ncells.append(6)
layer_mat_ids.append(1001)
z_curr -= 0.06

# Catotelm (0.27 m, 1 cm cells)
layer_types.append("constant")
layer_data.append(0.27)
layer_ncells.append(27)
layer_mat_ids.append(1002)
z_curr -= 0.27

# Mineral Layer Part 1 (0.6 m, 1 cm cells)
layer_types.append("constant")
layer_data.append(0.6)
layer_ncells.append(60)
layer_mat_ids.append(1003)
z_curr -= 0.6

# Mineral Layer Part 2: 50 cells with gradually increasing dz
dz = 0.02
for _ in range(50):
    dz *= 1.05
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1003)
    z_curr -= dz

# Mineral Layer Part 3: Continue exponential growth until dz reaches ~1.8 m
while dz < 1.8:
    dz *= 1.05
    layer_types.append("constant")
    layer_data.append(dz)
    layer_ncells.append(1)
    layer_mat_ids.append(1003)
    z_curr -= dz

# Mineral Layer Part 4: Constant 2.0 m cells until bottom of domain (40 m total)
remaining_depth = 40 + z_curr  # z_curr is negative
layer_types.append("constant")
layer_data.append(remaining_depth)
layer_ncells.append(int(round(remaining_depth / 2.0)))
layer_mat_ids.append(1003)

# ----------------------
# GENERATE 3D MESH AND WRITE
# ----------------------

m3 = meshing_ats.Mesh3D.extruded_Mesh2D(
    m2,
    layer_types,
    layer_data,
    layer_ncells,
    layer_mat_ids
)

output_file = "../data/mesh_transect_imnavait.exo"
m3.write_exodus(output_file)
print(f"Mesh written to: {output_file}")