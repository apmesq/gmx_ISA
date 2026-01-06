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
        print("Usage: python3 isa.py <Coulomb_energy>.xvg <LJ_energy>.xvg")
        sys.exit(1)

    coulomb_file = sys.argv[1]
    lj_file = sys.argv[2]

    if not os.path.isfile(coulomb_file):
        print(f"File not found: {coulomb_file}")
        sys.exit(1)

    if not os.path.isfile(lj_file):
        print(f"File not found: {lj_file}")
        sys.exit(1)

    # Read energies
    coulomb = read_xvg_energy(coulomb_file)
    lj = read_xvg_energy(lj_file)

    if len(coulomb) != len(lj):
        print("Error: Coulomb and LJ files must have the same number of points.")
        sys.exit(1)

    # Interaction energy
    interaction_energy = coulomb + lj
    N = len(interaction_energy)

    # Statistics
    mean_E = np.mean(interaction_energy)
    std_E = np.std(interaction_energy, ddof=1)

    if std_E == 0.0:
        print("Error: Standard deviation is zero. ISA is undefined.")
        sys.exit(1)

    # Errors
    mean_error = std_E / np.sqrt(N)
    std_error = std_E / np.sqrt(2 * (N - 1))

    ISA = mean_E / std_E

    ISA_error = abs(ISA) * np.sqrt(
        (mean_error / mean_E) ** 2 +
        (std_error / std_E) ** 2
    )

    # Output
    output_file = "ISA_results.dat"

    with open(output_file, "w") as f:
        f.write("Interaction Stability Analysis (ISA)\n\n")

        f.write("Definitions:\n")
        f.write("  E_int = E_Coul + E_LJ\n")
        f.write("  ISA   = <E_int> / sigma(E_int)\n\n")

        f.write(f"Number of frames             = {N}\n\n")

        f.write(f"Mean interaction energy      = {mean_E:.6f} +/- {mean_error:.6f}\n")
        f.write(f"Std interaction energy       = {std_E:.6f} +/- {std_error:.6f}\n")
        f.write(f"ISA                          = {ISA:.6f} +/- {ISA_error:.6f}\n")

    print("ISA analysis with uncertainties completed.")
    print(f"Results written to: {output_file}")


if __name__ == "__main__":
    main()

