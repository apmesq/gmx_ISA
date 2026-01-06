#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


# ============================================================
# Functions
# ============================================================
def read_xvg_energy(filename):
    """
    Reads a GROMACS .xvg file containing time vs energy.
    Returns only the energy values.
    """
    energy = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                energy.append(float(parts[1]))

    return np.array(energy)


def linear_regression(x, y):
    """
    Performs linear regression y = a*x + b.
    Returns slope, intercept, R^2, slope error, intercept error, fitted y.
    """
    x = x.reshape(-1, 1)
    model = LinearRegression()
    model.fit(x, y)

    y_fit = model.predict(x)

    slope = model.coef_[0]
    intercept = model.intercept_
    r2 = r2_score(y, y_fit)

    n = len(y)
    residuals = y - y_fit
    s_err = np.sqrt(np.sum(residuals**2) / (n - 2))

    x_var = np.sum((x.flatten() - np.mean(x))**2)

    slope_err = s_err / np.sqrt(x_var)
    intercept_err = s_err * np.sqrt(1.0 / n + (np.mean(x)**2) / x_var)

    return slope, intercept, r2, slope_err, intercept_err, y_fit


def plot_correlation(x, y, fit, outfile):
    """
    Plots Coulomb vs Lennard-Jones energy correlation.
    """
    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, s=15, alpha=0.6, label="Data")
    plt.plot(x, fit, linestyle="--", linewidth=2, label="Linear fit")
    plt.xlabel("Lennard-Jones Energy")
    plt.ylabel("Coulomb Energy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


# ============================================================
# Main
# ============================================================
def main():

    if len(sys.argv) != 3:
        print("Usage: python3 int.py <Coulomb_energy>.xvg <LJ_energy>.xvg")
        sys.exit(1)

    coulomb_file = sys.argv[1]
    lj_file = sys.argv[2]

    if not os.path.isfile(coulomb_file) or not os.path.isfile(lj_file):
        print("Error: One or more input files not found.")
        sys.exit(1)

    base = "interaction_energy"

    # Read energies
    coulomb = read_xvg_energy(coulomb_file)
    lj = read_xvg_energy(lj_file)

    if len(coulomb) != len(lj):
        print("Error: Coulomb and LJ files have different number of points.")
        sys.exit(1)

    # ========================================================
    # Statistics
    # ========================================================
    coulomb_mean = np.mean(coulomb)
    coulomb_std = np.std(coulomb, ddof=1)

    lj_mean = np.mean(lj)
    lj_std = np.std(lj, ddof=1)

    total = coulomb + lj
    total_mean = np.mean(total)
    total_std = np.std(total, ddof=1)

    # ========================================================
    # Coulomb vs LJ linear regression
    # ========================================================
    slope, intercept, r2, slope_err, intercept_err, fit = linear_regression(lj, coulomb)

    # ========================================================
    # Plot
    # ========================================================
    plot_correlation(lj, coulomb, fit, f"{base}_coulomb_vs_lj.png")

    # ========================================================
    # Output .dat file
    # ========================================================
    with open(f"{base}_analysis.dat", "w") as f:
        f.write("Interaction energy correlation analysis\n\n")

        f.write("Coulomb energy:\n")
        f.write(f"  Mean = {coulomb_mean:.6f}\n")
        f.write(f"  Std  = {coulomb_std:.6f}\n\n")

        f.write("Lennard-Jones energy:\n")
        f.write(f"  Mean = {lj_mean:.6f}\n")
        f.write(f"  Std  = {lj_std:.6f}\n\n")

        f.write("Total interaction energy (Coulomb + LJ):\n")
        f.write(f"  Mean = {total_mean:.6f}\n")
        f.write(f"  Std  = {total_std:.6f}\n\n")

        f.write("Coulomb vs Lennard-Jones linear regression:\n")
        f.write("  Model: E_Coulomb = a * E_LJ + b\n")
        f.write(f"  Slope (a)        = {slope:.6e} +/- {slope_err:.6e}\n")
        f.write(f"  Intercept (b)    = {intercept:.6e} +/- {intercept_err:.6e}\n")
        f.write(f"  R^2              = {r2:.6f}\n")

    print("Analysis completed successfully.")
    print("Generated files:")
    print(f"  {base}_analysis.dat")
    print(f"  {base}_coulomb_vs_lj.png")


if __name__ == "__main__":
    main()

