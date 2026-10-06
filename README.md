# Imnavait Creek rain-variability ensemble: model-data archive

Model inputs, outputs and analysis code supporting

> Mukherjee, N., Gao, B., Coon, E. T., Shuai, P., Hill, D., Neilson, B. T., Cory, R. M., Kling, G. W., Chen, J., & Cardenas, M. B.
> *The Effects of Precipitation Variability on Supra-Permafrost Thermal Hydrology: Flashy Hillslopes with More Outflow.*
> Submitted to *Water Resources Research*.

Simulations use the Advanced Terrestrial Simulator (ATS) v1.6-dev, commit `414f81aa` (https://github.com/amanzi/ats), on a 2D hillslope–riparian transect (109 m × 40 m) at Imnavait Creek, Alaska. The meteorological forcing is the 40-year (1981–2020) day-of-year mean of the NASA DAAC ABoVE SnowModel dataset (Liston et al., 2023). Summer rain is replaced by synthetic lognormal ensembles with nine levels of variability (σ1–σ9) and identical seasonal rain totals.

## Contents

| Path | Contents |
|---|---|
| `data/` | Meshes, forcing files, column initial conditions, and the Imnavait Creek weir discharge record |
| `inputs/` | ATS input files (`.xml`), one per run |
| `runs/` | Model output: final checkpoint and daily observations per run, plus trimmed visualization output for Figures 3–4 |
| `scripts/` | Mesh, forcing and initial-condition generators, the visualization trimming tool, and `ICVar_Figures.ipynb` |
| `flmd.csv` | File-level metadata (ESS-DIVE) |
| `dd.csv` | Data dictionary for the tabular files (ESS-DIVE) |

## σ labels

Run names, input files and forcing columns use a 0-based index; the manuscript uses σ1–σ9.

| Run index (`sigma{N}`, forcing column `rain precipitation N`) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|---|
| Manuscript | σ1 | σ2 | σ3 | σ4 | **σ5 (control)** | σ6 | σ7 | σ8 | σ9 |
| Lognormal shape parameter | 0.58 | 0.68 | 0.78 | 1.08 | 1.48 | 1.88 | 2.18 | 2.28 | 2.38 |

σ5 is the lognormal fit to observed summer rain (μ = 0.40, σ = 1.48); the others are offset by ±0.4, ±0.7, ±0.8 and ±0.9 with the mean held fixed.

## Simulation workflow

| Step | Run | Input | Initial condition | Rain forcing | Period |
|---|---|---|---|---|---|
| 1 | `run01_column_freezeup` | `column_freezeup_imnavait.xml` | see input file | — | freeze-up |
| 2 | `run02_column_spinup_smoothedmet` | `column_spinup_imnavait_smoothedMet.xml` | `data/column_data_freezeup.h5` | `data/forcing_smoothed4spinUp.h5` | spinup |
| 3 | `run03_transect_base_sigma{N}_6yr` | `transect_base_sigma{N}_6yr.xml` | `data/column_data_spinup_smoothedMet.h5` | `synthetic_rainfall_ensemble_6yr.h5` | 6 yr (days 1–2190) |
| 4 | `run04_transect_extra_sigma{N}_3yr` | `transect_extra_sigma{N}_3yr.xml` | restart `run03_…/checkpoint_final.h5` | `synthetic_rainfall_ensemble_9yr.h5` (years 7–9) | 3 yr (days 2190–3285) |
| 5 | `run05_transect_seeds_sigma{N}_4yr` | `transect_seeds_sigma{N}_4yr.xml` | restart `run03_…/checkpoint_final.h5` | `synthetic_rainfall_ensemble_4yr_new.h5` | 4 yr (days 0–1460) |

Each transect run is started from its own directory under `runs/` with `ats --xml_file=../../inputs/<file>.xml`; all paths in the input files are relative to that directory.

### Ensemble realizations

Each σ level yields 12 one-year thaw-season realizations:

| Realizations | Source | Random seeds |
|---|---|---|
| a–e | `run03` years 2–6 (year 1 is warm-up) | shared across σ levels |
| f–h | `run04` years 1–3 | shared across σ levels |
| i–l | `run05` years 1–4 | unique to each σ level |

### Preprocessing scripts (`scripts/`)

| Script | Produces |
|---|---|
| `gen_column_imnavait.py` | `data/mesh_column_imnavait.exo` |
| `gen_imnavait_transect.py` | `data/mesh_transect_imnavait.exo` |
| `smooth_met.py` | Smoothed single-year spinup forcing (`data/forcing_smoothed4spinUp.h5`) |
| `column_data.py` | `column_data.h5` from a finished column run (pressure/temperature vs. depth). Run in `run01` → `data/column_data_freezeup.h5`; run in `run02` → `data/column_data_spinup_smoothedMet.h5` |
| `gen_ensemble_rain_forcings.py` | The three `synthetic_rainfall_ensemble_*.h5` files from `forcing_daac_ABoVE_snowmodel_1981_2020.h5` and the weir discharge record (bit-for-bit reproduction of the archived files) |
| `trim_vis_output.py` | The trimmed `ats_vis_*.h5` files in `run03_transect_base_sigma{0,4,8}_6yr` (output frames 1293–1295 from the full ~10 GB visualization output) |

The mesh and initial-condition scripts require the ATS Python tools (`$ATS_SRC_DIR/tools/utils`: `ats_xdmf`, `meshing_ats`) and SEACAS `exodus`.

## Outputs (`runs/`)

* `checkpoint_final.h5`: ATS restart file at the end of the run
* `observation_2D.dat` (transect) / `observation_columnSpinup.dat` (column): daily domain-integrated or averaged observations, comma-separated, with `#` header lines describing each column's region and reduction. Column definitions are in `dd.csv`. Convert `[mol]` water quantities to m³ with 55 500 mol m⁻³.
* `ats_vis_data.h5`, `ats_vis_mesh.h5`, `ats_vis_surface_data.h5`, `ats_vis_surface_mesh.h5` (σ1, σ5, σ9 base runs only): subsurface saturation, temperature, pressure and Darcy flux, plus surface elevation, ponded depth and rain, at simulation days 1294–1296 (18–20 July of year 4). Readable with `ats_xdmf.VisFile(<run dir>)`.

The full visualization output (2190 frames per base run) is not archived.

## Figures

`scripts/ICVar_Figures.ipynb` regenerates the data panels of main-text Figures 1–8 from this package into `scripts/figures/` (≈2 min). Figure 9 is conceptual. Python ≥ 3.10 with numpy, pandas, scipy, matplotlib, seaborn and h5py, plus `ats_xdmf` from the ATS source tree (located through `$ATS_SRC_DIR`).

## External data

* **Meteorological forcing**: `data/forcing_daac_ABoVE_snowmodel_1981_2020.h5`, daily SnowModel output at Imnavait Creek from Liston, G. E., et al. (2023), NASA ORNL DAAC ABoVE SnowModel dataset, converted to ATS forcing format.
* **Discharge**: `data/ImnavaitCr_Historical_Discharge_1985_2017_2018Aug6.csv`, from Kane, D. L., Hinzman, L. D., Stuefer, S. L., et al., *Hydrographic Data, Imnavait Creek Watershed, Alaska, 1985–2017*, Arctic Data Center. Used only to define the mean summer window (freshet to end of streamflow).
