#!/usr/bin/env python3

import sys
import os
import numpy as np


def read_xvg_energy(filename):
    """
    Reads a GROMACS .xvg file containing energy vs time.
    Returns a NumPy array with the energy values.
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


def main():

    if len(sys.argv) != 3:
        print("Usage: python3 int.py <Coulomb_energy>.xvg <LJ_energy>.xvg")
        sys.exit(1)

    coulomb_file = sys.argv[1]
    lj_file = sys.argv[2]

    if not os.path.isfile(coulomb_file):
        print(f"File not found: {coulomb_file}")
        sys.exit(1)

    if not os.path.isfile(lj_file):
        print(f"File not found: {lj_file}")
        sys.exit(1)

    # Read energy data
    coulomb = read_xvg_energy(coulomb_file)
    lj = read_xvg_energy(lj_file)

    if len(coulomb) != len(lj):
        print("Error: Coulomb and LJ files do not have the same number of data points.")
        sys.exit(1)

    # Total interaction energy
    total = coulomb + lj

    # Statistics
    coulomb_mean = np.mean(coulomb)
    coulomb_std = np.std(coulomb, ddof=1)

    lj_mean = np.mean(lj)
    lj_std = np.std(lj, ddof=1)

    total_mean = np.mean(total)
    total_std = np.std(total, ddof=1)

    # Output file
    output_file = "interaction_energy_statistics.dat"

    with open(output_file, "w") as f:
        f.write("Interaction energy statistics (from GROMACS XVG files)\n\n")

        f.write("Coulomb energy:\n")
        f.write(f"  Mean = {coulomb_mean:.6f}\n")
        f.write(f"  Std  = {coulomb_std:.6f}\n\n")

        f.write("Lennard-Jones energy:\n")
        f.write(f"  Mean = {lj_mean:.6f}\n")
        f.write(f"  Std  = {lj_std:.6f}\n\n")

        f.write("Total interaction energy (Coulomb + LJ):\n")
        f.write(f"  Mean = {total_mean:.6f}\n")
        f.write(f"  Std  = {total_std:.6f}\n")

    print("Analysis completed successfully.")
    print(f"Results written to: {output_file}")


if __name__ == "__main__":
    main()

