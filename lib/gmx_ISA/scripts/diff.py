#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import t


# ============================================================
# Constants
# ============================================================
DIMENSIONALITY_FACTOR = 6.0  # 3D diffusion: MSD = 6Dt


# ============================================================
# Functions
# ============================================================
def read_xvg(filename):
    """
    Reads a GROMACS .xvg file and returns x and y arrays.
    Ignores comment lines starting with @ or #.
    """
    x, y = [], []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))

    return np.array(x), np.array(y)


def linear_regression(x, y):
    """
    Performs linear regression y = ax + b.
    Returns slope, intercept, R^2, and standard error of slope.
    """
    x = x.reshape(-1, 1)
    model = LinearRegression()
    model.fit(x, y)
    y_pred = model.predict(x)

    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = r2_score(y, y_pred)

    # Standard error of the slope
    n = len(y)
    residuals = y - y_pred
    s_err = np.sqrt(np.sum(residuals**2) / (n - 2))
    x_var = np.sum((x.flatten() - np.mean(x))**2)
    slope_err = s_err / np.sqrt(x_var)

    return slope, intercept, r2, slope_err


def find_best_linear_region(x, y, min_fraction=0.5):
    """
    Finds the most linear continuous interval containing at least
    min_fraction of the points, maximizing R^2.
    """
    n = len(x)
    min_points = int(np.ceil(min_fraction * n))

    best_r2 = -np.inf
    best = None

    for start in range(0, n - min_points + 1):
        for end in range(start + min_points, n + 1):
            x_win = x[start:end]
            y_win = y[start:end]

            slope, intercept, r2, slope_err = linear_regression(x_win, y_win)

            if r2 > best_r2:
                best_r2 = r2
                best = {
                    "start": start,
                    "end": end,
                    "slope": slope,
                    "intercept": intercept,
                    "r2": r2,
                    "slope_err": slope_err
                }

    return best


# ============================================================
# Main
# ============================================================
def main():

    if len(sys.argv) != 2:
        print("Usage: python3 diff.py <msd_file>.xvg")
        sys.exit(1)

    filename = sys.argv[1]

    if not os.path.isfile(filename):
        print(f"File not found: {filename}")
        sys.exit(1)

# PEGA APENAS O NOME DO ARQUIVO, IGNORANDO O DIRETÓRIO ANTERIOR
    base = os.path.splitext(os.path.basename(filename))[0]

    # Read MSD data
    time, msd = read_xvg(filename)

    # ========================================================
    # Full data regression
    # ========================================================
    slope_all, intercept_all, r2_all, slope_err_all = linear_regression(time, msd)

    D_all = slope_all / DIMENSIONALITY_FACTOR
    D_all_err = slope_err_all / DIMENSIONALITY_FACTOR

    # ========================================================
    # Best continuous region (>= 50%)
    # ========================================================
    best = find_best_linear_region(time, msd, min_fraction=0.5)

    s, e = best["start"], best["end"]

    time_best = time[s:e]
    msd_best = msd[s:e]

    D_best = best["slope"] / DIMENSIONALITY_FACTOR
    D_best_err = best["slope_err"] / DIMENSIONALITY_FACTOR

    # ========================================================
    # Plot: full MSD
    # ========================================================
    plt.figure(figsize=(7, 5))
    plt.plot(time, msd, label="MSD (all points)")
    plt.plot(time, slope_all * time + intercept_all,
             linestyle="--", label="Linear fit (all points)")
    plt.xlabel("Time")
    plt.ylabel("MSD")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{base}_msd_all.png", dpi=300)
    plt.close()

    # ========================================================
    # Plot: best linear region
    # ========================================================
    plt.figure(figsize=(7, 5))
    plt.plot(time, msd, color="gray", alpha=0.4, label="MSD (all points)")
    plt.plot(time_best, msd_best, label="Best linear region")
    plt.plot(time_best,
             best["slope"] * time_best + best["intercept"],
             linestyle="--", label="Linear fit (best region)")
    plt.xlabel("Time")
    plt.ylabel("MSD")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{base}_msd_best_region.png", dpi=300)
    plt.close()

    # ========================================================
    # Write results file (.dat)
    # ========================================================
    with open(f"{base}_diffusion_results.dat", "w") as f:
        f.write("Diffusion coefficient analysis from MSD\n\n")

        f.write("Full data regression:\n")
        f.write(f"  Slope              = {slope_all:.6e}\n")
        f.write(f"  Diffusion (D)      = {D_all:.6e} +/- {D_all_err:.6e}\n")
        f.write(f"  R^2                = {r2_all:.6f}\n\n")

        f.write("Best continuous region (>= 50% of points):\n")
        f.write(f"  Index range        = {s} to {e-1}\n")
        f.write(f"  Time range         = {time[s]:.6f} to {time[e-1]:.6f}\n")
        f.write(f"  Slope              = {best['slope']:.6e}\n")
        f.write(f"  Diffusion (D)      = {D_best:.6e} +/- {D_best_err:.6e}\n")
        f.write(f"  R^2                = {best['r2']:.6f}\n")

    # ========================================================
    # Export best region to CSV
    # ========================================================
    df_best = pd.DataFrame({
        "time": time_best,
        "msd": msd_best
    })

    df_best.to_csv(f"{base}_best_linear_region.csv", index=False)

    print("Analysis completed successfully.")
    print(f"Results written to:")
    print(f"  {base}_diffusion_results.dat")
    print(f"  {base}_best_linear_region.csv")
    print(f"  {base}_msd_all.png")
    print(f"  {base}_msd_best_region.png")


if __name__ == "__main__":
    main()

