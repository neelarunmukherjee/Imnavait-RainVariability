"""
trim_vis_output.py

Extract a small subset of ATS visualization output (ats_vis_*.h5) for archiving.

The full ats_vis_data.h5 of each 2D transect run is ~10 GB (2190 output cycles x
29 fields). Figure 3 of the manuscript only needs the subsurface saturation and
Darcy velocity fields at a few output frames (around a summer storm event), plus
the static surface elevation. This script copies just those variables and frames
into new HDF5 files with the same layout ATS writes, so they can be read with
ats_xdmf.VisFile exactly like the original output.

Usage:
  python trim_vis_output.py SRC_RUN_DIR DST_RUN_DIR [--frames 1293 1294 1295]

  SRC_RUN_DIR : original run directory containing ats_vis_*.h5
  DST_RUN_DIR : archive run directory (receives ats_vis_*.h5 subsets)
  --frames    : 0-based indices into the original (sorted) output cycle list

Notes:
  - Frame indices refer to the ORIGINAL output; in the trimmed files the same
    frames are indices 0..len(frames)-1. Each cycle keeps its 'Time' attribute,
    so select frames by time when reading the trimmed files.
  - Mesh files (ats_vis_mesh.h5, ats_vis_surface_mesh.h5) are copied with the
    mesh group re-keyed to the first trimmed cycle (the mesh is static).
"""

import argparse
import os

import h5py

SUBSURFACE_VARS = [
    "cell_volume",          # required by ats_xdmf.VisFile.loadMesh
    "saturation_liquid",
    "saturation_ice",
    "temperature",
    "pressure",
    "darcy_velocity.0",
    "darcy_velocity.1",
    "darcy_velocity.2",
]
SURFACE_VARS = [
    "surface-cell_volume",  # required by ats_xdmf.VisFile.loadMesh
    "surface-elevation",
    "surface-ponded_depth",
    "surface-precipitation_rain",
]


def trim(src_file, dst_file, variables, cycles):
    """Copy `variables` at `cycles` from src_file to dst_file, keeping attributes."""
    with h5py.File(src_file, "r") as src, h5py.File(dst_file, "w") as dst:
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        for var in variables:
            if var not in src:
                raise KeyError(f"{var} not found in {src_file}")
            g = dst.create_group(var)
            for k, v in src[var].attrs.items():
                g.attrs[k] = v
            for c in cycles:
                src.copy(src[var][c], g, name=c)


def copy_mesh(src_file, dst_file, first_cycle):
    """Copy the (static) mesh, keyed by `first_cycle`.

    ats_xdmf.VisFile looks up the mesh under the first data cycle, so the mesh
    group written at the original first cycle is stored under the first
    trimmed cycle instead.
    """
    with h5py.File(src_file, "r") as src, h5py.File(dst_file, "w") as dst:
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        src.copy(src[list(src.keys())[0]], dst, name=first_cycle)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("src_run_dir")
    p.add_argument("dst_run_dir")
    p.add_argument("--frames", type=int, nargs="+", default=[1293, 1294, 1295])
    args = p.parse_args()

    os.makedirs(args.dst_run_dir, exist_ok=True)
    src_data = os.path.join(args.src_run_dir, "ats_vis_data.h5")

    with h5py.File(src_data, "r") as f:
        all_cycles = sorted(f[SUBSURFACE_VARS[0]].keys(), key=int)
        cycles = [all_cycles[i] for i in args.frames]
        times = [f[SUBSURFACE_VARS[0]][c].attrs["Time"] for c in cycles]

    trim(src_data, os.path.join(args.dst_run_dir, "ats_vis_data.h5"),
         SUBSURFACE_VARS, cycles)
    trim(os.path.join(args.src_run_dir, "ats_vis_surface_data.h5"),
         os.path.join(args.dst_run_dir, "ats_vis_surface_data.h5"),
         SURFACE_VARS, cycles)
    for mesh in ["ats_vis_mesh.h5", "ats_vis_surface_mesh.h5"]:
        copy_mesh(os.path.join(args.src_run_dir, mesh),
                  os.path.join(args.dst_run_dir, mesh), cycles[0])

    for i, c, t in zip(args.frames, cycles, times):
        print(f"frame {i}: cycle {c}, time {t}")


if __name__ == "__main__":
    main()
