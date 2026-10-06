"""
gen_ensemble_rain_forcings.py

Generate synthetic ensemble rainfall forcing files for ATS transect simulations.

Reproduces the logic from three notebooks:
  - make_forcings_vcariability_DAAC_final_emsemble_rain_0.ipynb  -> 6yr ensemble
  - make_forcings_vcariability_DAAC_final_emsemble_rain.ipynb    -> 9yr ensemble
  - make_forcings_vcariability_DAAC_final_emsemble_rain_seeds.ipynb -> 4yr_new

Outputs (written to OUTPUT_DIR):
  - synthetic_rainfall_ensemble_6yr.h5         (9 sigmas x 6 reps, a-f)
  - synthetic_rainfall_ensemble_9yr.h5         (9 sigmas x 9 reps, a-i)
  - synthetic_rainfall_ensemble_4yr_new.h5     (9 sigmas x 4 reps, j-m)

Usage:
  python gen_ensemble_rain_forcings.py
  python gen_ensemble_rain_forcings.py --output-dir /path/to/data

Key parameters:
  - Forcing data: DAAC ABoVE snowmodel daily met (1981-2020)
  - Discharge data: Imnavait Creek historical discharge (1985-2017)
  - Summer season defined by mean freshet DOY and end-of-discharge DOY
  - Lognormal fit to combined summer rain; mean held fixed, sigma varied
  - Conversion: 1.5741e-8 (mm/d -> m/s)
"""

import argparse
import os
import numpy as np
import pandas as pd
import h5py
from scipy.optimize import curve_fit


# ── Paths ────────────────────────────────────────────────────────────────────

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "../data")

FORCING_FILE = os.path.join(DATA_DIR, "forcing_daac_ABoVE_snowmodel_1981_2020.h5")
DISCHARGE_CSV = os.path.join(DATA_DIR, "ImnavaitCr_Historical_Discharge_1985_2017_2018Aug6.csv")


# ── Constants ────────────────────────────────────────────────────────────────

CONVERSION_FACTOR = 1.5741e-8   # mm/d -> m/s


# ── Helper functions ─────────────────────────────────────────────────────────

def lognormal_pdf(x, mu, sigma):
    """Lognormal probability density function."""
    return (1.0 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(
        -((np.log(x) - mu) ** 2) / (2 * sigma ** 2)
    )


def generate_fixed_mean_precip(sigma, desired_mean, days):
    """
    Draw lognormal samples then rescale to exactly desired_mean.
    Caller must set np.random.seed() before calling for reproducibility.
    """
    mu = np.log(desired_mean) - (sigma ** 2) / 2
    precip = np.random.lognormal(mean=mu, sigma=sigma, size=days)
    precip *= desired_mean / np.mean(precip)
    return precip


# ── Load DAAC forcing ────────────────────────────────────────────────────────

def load_daac_forcing(forcing_file):
    """Load DAAC ABoVE snowmodel HDF5 into a dict of 1-D numpy arrays."""
    with h5py.File(forcing_file, "r") as f:
        data = {key: f[key][:] for key in f.keys()}
    return data


# ── Main ─────────────────────────────────────────────────────────────────────

def main(output_dir):
    os.makedirs(output_dir, exist_ok=True)

    # ── 1. Load forcing ──────────────────────────────────────────────────────
    print("Loading DAAC forcing...")
    raw = load_daac_forcing(FORCING_FILE)

    dates = pd.date_range(start="1981-01-01", end="2020-08-31", freq="D")

    df_forcing = pd.DataFrame({
        "Precipitation Rain [m s^-1]":          raw["Precipitation Rain [m s^-1]"],
        "Precipitation Snow [m SWE s^-1]":      raw["Precipitation Snow [m SWE s^-1]"],
        "air temperature [K]":                  raw["air temperature [K]"],
        "incoming shortwave radiation [W m^-2]": raw["incoming shortwave radiation [W m^-2]"],
        "vapor pressure air [Pa]":              raw["vapor pressure air [Pa]"],
        "wind speed [m s^-1]":                  raw["wind speed [m s^-1]"],
    })
    df_forcing["Date"] = dates
    df_forcing["Precipitation Rain [mm/d]"] = df_forcing["Precipitation Rain [m s^-1]"] * 8.64e7
    df_forcing["air temperature [C]"] = df_forcing["air temperature [K]"] - 273.15
    df_forcing["Year"] = df_forcing["Date"].dt.year
    df_forcing["DOY"] = df_forcing["Date"].dt.day_of_year
    years = df_forcing["Year"].unique()

    # ── 2. DOY climatological means ──────────────────────────────────────────
    print("Computing DOY climatological means...")
    grouped = df_forcing.groupby("DOY")
    means = {
        "rain precipitation [m/s]":              grouped["Precipitation Rain [m s^-1]"].mean(),
        "snow precipitation [m SWE/s]":          grouped["Precipitation Snow [m SWE s^-1]"].mean(),
        "air temperature [K]":                   grouped["air temperature [K]"].mean(),
        "incoming shortwave radiation [W m^-2]": grouped["incoming shortwave radiation [W m^-2]"].mean(),
        "vapor pressure air [Pa]":               grouped["vapor pressure air [Pa]"].mean(),
        "wind speed [m/s]":                      grouped["wind speed [m s^-1]"].mean(),
    }
    means_df = pd.DataFrame(means)
    means_df = means_df.iloc[:-1]          # drop DOY 366 (leap-year artefact)
    means_df.index = range(len(means_df))  # reset to 0-364

    # ── 3. Determine summer season DOYs from discharge ───────────────────────
    print("Computing summer season from discharge record...")
    df_dis = pd.read_csv(DISCHARGE_CSV)
    df_clean = df_dis[["Date", "Discharge"]].copy()
    df_clean["Date"] = pd.to_datetime(df_clean["Date"], errors="coerce")
    df_clean["Discharge"] = pd.to_numeric(df_clean["Discharge"], errors="coerce")
    df_clean.loc[df_clean["Discharge"] == 6999, "Discharge"] = np.nan
    df_clean.dropna(subset=["Date", "Discharge"], inplace=True)
    df_clean["Year"] = df_clean["Date"].dt.year
    df_clean["Month"] = df_clean["Date"].dt.month

    freshet_dates = []
    for year, group in df_clean.groupby("Year"):
        may_data = group[group["Month"] == 5]
        summer_data = group[(group["Month"] >= 5) & (group["Month"] <= 9)]
        if not may_data.empty:
            freshet_date = may_data.loc[may_data["Discharge"].idxmax(), "Date"]
        else:
            freshet_date = np.nan
        summer_nonzero = summer_data[summer_data["Discharge"] > 0]
        end_of_discharge_date = summer_nonzero["Date"].iloc[-1] if not summer_nonzero.empty else np.nan
        freshet_dates.append(
            {"Year": year, "Freshet Date": freshet_date, "End of Discharge": end_of_discharge_date}
        )

    result_df = pd.DataFrame(freshet_dates)
    mean_freshet_doy = result_df["Freshet Date"].dropna().dt.dayofyear.mean()
    mean_end_doy = result_df["End of Discharge"].dropna().dt.dayofyear.mean()
    mean_freshet_date = pd.to_datetime("2020-01-01") + pd.to_timedelta(mean_freshet_doy - 1, unit="D")
    mean_end_date = pd.to_datetime("2020-01-01") + pd.to_timedelta(mean_end_doy - 1, unit="D")
    start_doy = int(mean_freshet_date.day_of_year)
    end_doy = int(mean_end_date.day_of_year)
    days = end_doy - start_doy + 1
    print(f"  start_doy={start_doy}, end_doy={end_doy}, days={days}")

    # ── 4. Fit lognormal to combined summer rain ─────────────────────────────
    print("Fitting lognormal to summer precipitation...")
    rain_col = "Precipitation Rain [mm/d]"
    combined_data = df_forcing[(df_forcing["DOY"] >= start_doy) & (df_forcing["DOY"] <= end_doy)]

    hist_vals, bin_edges = np.histogram(combined_data[rain_col], bins=50, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    valid = bin_centers > 0
    bin_centers = bin_centers[valid]
    hist_vals = hist_vals[valid]

    rain_pos = combined_data[rain_col][combined_data[rain_col] > 0]
    initial_guess = [np.mean(np.log(rain_pos)), np.std(np.log(rain_pos))]
    popt, _ = curve_fit(lognormal_pdf, bin_centers, hist_vals, p0=initial_guess)
    fitted_mu, fitted_sigma = popt

    logn_mean = np.exp(fitted_mu + (fitted_sigma ** 2) / 2)
    mean_summer_precip = np.mean(combined_data[rain_col])
    print(f"  fitted mu={fitted_mu:.4f}, sigma={fitted_sigma:.4f}")
    print(f"  lognormal mean={logn_mean:.4f} mm/d, sample mean={mean_summer_precip:.4f} mm/d")

    # ── 5. Define 9 sigma levels ─────────────────────────────────────────────
    s = fitted_sigma
    all_sigmas = sorted([
        s - 0.9, s - 0.8, s - 0.7, s - 0.4,
        s,
        s + 0.4, s + 0.7, s + 0.8, s + 0.9,
    ])
    print(f"  sigma levels: {[f'{v:.4f}' for v in all_sigmas]}")

    # ── 6. Generate synthetic precip series ──────────────────────────────────
    print("Generating synthetic precipitation (reps a-i)...")
    date_idx_summer = pd.date_range(start="2094-05-20", periods=days)

    # Reps a-i: seeds [50,46,47,48,49,45,44,51,52]
    seeds_orig = [50, 46, 47, 48, 49, 45, 44, 51, 52]
    reps_orig = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]
    synthetic_precip = {}
    for sigma_val, sigma in enumerate(all_sigmas):
        for seed, rep in zip(seeds_orig, reps_orig):
            np.random.seed(seed)
            synthetic_precip[f"sig{sigma_val}{rep}"] = generate_fixed_mean_precip(
                sigma, mean_summer_precip, days
            )
    df_synthetic_precip = pd.DataFrame(synthetic_precip, index=date_idx_summer)

    # Reps j-m: new unique seeds range(53, 89) cycling through sigma_val × rep
    print("Generating synthetic precipitation (reps j-m)...")
    new_seeds = list(range(53, 89))  # 36 unique seeds (9 sigmas × 4 reps)
    reps_new = ["j", "k", "l", "m"]
    synthetic_precip_new = {}
    seed_idx = 0
    for sigma_val, sigma in enumerate(all_sigmas):
        for rep in reps_new:
            np.random.seed(new_seeds[seed_idx])
            synthetic_precip_new[f"sig{sigma_val}{rep}"] = generate_fixed_mean_precip(
                sigma, mean_summer_precip, days
            )
            seed_idx += 1
    df_synthetic_precip_new = pd.DataFrame(synthetic_precip_new, index=date_idx_summer)

    # Re-index from start_doy so .loc[start_doy:end_doy] works
    df_synthetic_precip.index = range(start_doy, start_doy + days)
    df_synthetic_precip_new.index = range(start_doy, start_doy + days)

    # ── 7. DOY mean baseline forcing (1 year, 365 days) ──────────────────────
    print("Building caseA (climatological mean forcing)...")
    caseA_df = pd.DataFrame(index=np.arange(0, 365))
    caseA_df["air temperature [K]"] = means_df["air temperature [K]"].values
    caseA_df["vapor pressure air [Pa]"] = means_df["vapor pressure air [Pa]"].values
    caseA_df["incoming shortwave radiation [W m^-2]"] = means_df["incoming shortwave radiation [W m^-2]"].values
    caseA_df["wind speed [m/s]"] = means_df["wind speed [m/s]"].values
    caseA_df["snow precipitation [m SWE/s]"] = means_df["snow precipitation [m SWE/s]"].values
    caseA_df.loc[start_doy:end_doy, "snow precipitation [m SWE/s]"] = 0.0  # zero snow in summer

    # ── 8. Build case_ensemble helpers ───────────────────────────────────────

    def build_case_ensemble(df_precip, rep_list, n_years):
        """Build n_years-long daily forcing dataframe; inject summer rain per replicate."""
        ce = pd.DataFrame(index=np.arange(0, 365 * n_years))
        for sigma_val in range(9):
            col = f"rain precipitation {sigma_val} [m/s]"
            ce[col] = 0.0
            for rep_idx, rep in enumerate(rep_list):
                s_idx = start_doy + rep_idx * 365
                e_idx = end_doy + rep_idx * 365
                label = f"sig{sigma_val}{rep}"
                ce.loc[s_idx:e_idx, col] = (df_precip[label] * CONVERSION_FACTOR).to_numpy()
        return ce

    # ── 9. Write 6yr ensemble ────────────────────────────────────────────────
    print("Writing synthetic_rainfall_ensemble_6yr.h5 ...")
    case_ensemble_6yr = build_case_ensemble(df_synthetic_precip, reps_orig[:6], 6)
    case_means_6yr = pd.concat([caseA_df] * 6, ignore_index=True)

    dout_6yr = {}
    dout_6yr["time [s]"] = np.arange(0, len(case_ensemble_6yr)) * 86400
    dout_6yr["air temperature [K]"] = case_means_6yr["air temperature [K]"].astype(np.float32)
    dout_6yr["incoming shortwave radiation [W m^-2]"] = case_means_6yr["incoming shortwave radiation [W m^-2]"].astype(np.float32)
    dout_6yr["vapor pressure air [Pa]"] = case_means_6yr["vapor pressure air [Pa]"].astype(np.float32)
    dout_6yr["wind speed [m s^-1]"] = case_means_6yr["wind speed [m/s]"].astype(np.float32)
    dout_6yr["snow precipitation [m SWE s^-1]"] = case_means_6yr["snow precipitation [m SWE/s]"].astype(np.float32)
    for sigma_val in range(9):
        key_in = f"rain precipitation {sigma_val} [m/s]"
        key_out = f"rain precipitation {sigma_val} [m s^-1]"
        dout_6yr[key_out] = case_ensemble_6yr[key_in].astype(np.float32)

    out_path = os.path.join(output_dir, "synthetic_rainfall_ensemble_6yr.h5")
    with h5py.File(out_path, "w") as fid:
        for key, val in dout_6yr.items():
            fid.create_dataset(key, data=np.array(val, dtype=np.float32))
    print(f"  -> {out_path}  ({len(case_ensemble_6yr)} timesteps)")

    # ── 10. Write 9yr ensemble ───────────────────────────────────────────────
    print("Writing synthetic_rainfall_ensemble_9yr.h5 ...")
    case_ensemble_9yr = build_case_ensemble(df_synthetic_precip, reps_orig, 9)
    case_means_9yr = pd.concat([caseA_df] * 9, ignore_index=True)

    dout_9yr = {}
    dout_9yr["time [s]"] = np.arange(0, len(case_ensemble_9yr)) * 86400
    dout_9yr["air temperature [K]"] = case_means_9yr["air temperature [K]"].astype(np.float32)
    dout_9yr["incoming shortwave radiation [W m^-2]"] = case_means_9yr["incoming shortwave radiation [W m^-2]"].astype(np.float32)
    dout_9yr["vapor pressure air [Pa]"] = case_means_9yr["vapor pressure air [Pa]"].astype(np.float32)
    dout_9yr["wind speed [m s^-1]"] = case_means_9yr["wind speed [m/s]"].astype(np.float32)
    dout_9yr["snow precipitation [m SWE s^-1]"] = case_means_9yr["snow precipitation [m SWE/s]"].astype(np.float32)
    for sigma_val in range(9):
        key_in = f"rain precipitation {sigma_val} [m/s]"
        key_out = f"rain precipitation {sigma_val} [m s^-1]"
        dout_9yr[key_out] = case_ensemble_9yr[key_in].astype(np.float32)

    out_path = os.path.join(output_dir, "synthetic_rainfall_ensemble_9yr.h5")
    with h5py.File(out_path, "w") as fid:
        for key, val in dout_9yr.items():
            fid.create_dataset(key, data=np.array(val, dtype=np.float32))
    print(f"  -> {out_path}  ({len(case_ensemble_9yr)} timesteps)")

    # ── 11. Write 4yr_new ensemble ───────────────────────────────────────────
    print("Writing synthetic_rainfall_ensemble_4yr_new.h5 ...")
    case_ensemble_4yr = build_case_ensemble(df_synthetic_precip_new, reps_new, 4)
    case_means_4yr = pd.concat([caseA_df] * 4, ignore_index=True)

    dout_4yr = {}
    dout_4yr["time [s]"] = np.arange(0, len(case_ensemble_4yr)) * 86400
    dout_4yr["air temperature [K]"] = case_means_4yr["air temperature [K]"].astype(np.float32)
    dout_4yr["incoming shortwave radiation [W m^-2]"] = case_means_4yr["incoming shortwave radiation [W m^-2]"].astype(np.float32)
    dout_4yr["vapor pressure air [Pa]"] = case_means_4yr["vapor pressure air [Pa]"].astype(np.float32)
    dout_4yr["wind speed [m s^-1]"] = case_means_4yr["wind speed [m/s]"].astype(np.float32)
    dout_4yr["snow precipitation [m SWE s^-1]"] = case_means_4yr["snow precipitation [m SWE/s]"].astype(np.float32)
    for sigma_val in range(9):
        key_in = f"rain precipitation {sigma_val} [m/s]"
        key_out = f"rain precipitation {sigma_val} [m s^-1]"
        dout_4yr[key_out] = case_ensemble_4yr[key_in].astype(np.float32)

    out_path = os.path.join(output_dir, "synthetic_rainfall_ensemble_4yr_new.h5")
    with h5py.File(out_path, "w") as fid:
        for key, val in dout_4yr.items():
            fid.create_dataset(key, data=np.array(val, dtype=np.float32))
    print(f"  -> {out_path}  ({len(case_ensemble_4yr)} timesteps)")

    print("\nDone. All ensemble forcing files written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--output-dir",
        default=DATA_DIR,
        help="Directory to write output HDF5 files (default: ../data relative to script)",
    )
    args = parser.parse_args()
    main(os.path.abspath(args.output_dir))
