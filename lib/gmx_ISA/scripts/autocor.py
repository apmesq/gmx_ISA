#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def read_xvg(filename):
    """
    Reads a GROMACS .xvg file.
    Returns time array and data array.
    """
    time = []
    data = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                time.append(float(parts[0]))
                data.append(float(parts[1]))

    return np.array(time), np.array(data)


def autocorrelation_fft(x):
    """
    Computes normalized autocorrelation using FFT.
    """
    x = x - np.mean(x)
    n = len(x)

    fft_x = np.fft.fft(x, n=2*n)
    acf = np.fft.ifft(fft_x * np.conjugate(fft_x)).real
    acf = acf[:n]

    acf /= acf[0]
    return acf


def correlation_time(acf, dt):
    """
    Computes correlation time by integrating ACF
    until first zero crossing.
    """
    positive = np.where(acf > 0)[0]
    if len(positive) == 0:
        return 0.0

    k_cut = positive[-1]
    return np.sum(acf[:k_cut + 1]) * dt


def block_correlation_time(signal, time, n_blocks=5):
    """
    Estimates correlation time error using block averaging.
    """
    block_size = len(signal) // n_blocks
    tau_blocks = []

    for i in range(n_blocks):
        start = i * block_size
        end = (i + 1) * block_size

        block = signal[start:end]
        dt = np.mean(np.diff(time[start:end]))

        acf_block = autocorrelation_fft(block)
        tau_blocks.append(correlation_time(acf_block, dt))

    tau_blocks = np.array(tau_blocks)
    tau_mean = np.mean(tau_blocks)
    tau_err = np.std(tau_blocks, ddof=1) / np.sqrt(n_blocks)

    return tau_mean, tau_err, tau_blocks


def main():

    if len(sys.argv) != 3:
        print("Usage: python3 interaction_autocorr.py <Coulomb>.xvg <LJ>.xvg")
        sys.exit(1)

    coulomb_file = sys.argv[1]
    lj_file = sys.argv[2]

    if not os.path.isfile(coulomb_file) or not os.path.isfile(lj_file):
        print("Error: Input files not found.")
        sys.exit(1)

    # Read data
    time_c, coulomb = read_xvg(coulomb_file)
    time_lj, lj = read_xvg(lj_file)

    if len(coulomb) != len(lj):
        print("Error: Files must have the same length.")
        sys.exit(1)

    if not np.allclose(time_c, time_lj):
        print("Error: Time columns do not match.")
        sys.exit(1)

    time = time_c
    dt = np.mean(np.diff(time))

    # Interaction energy
    interaction_energy = coulomb + lj

    # Autocorrelation
    acf = autocorrelation_fft(interaction_energy)
    lags = np.arange(len(acf)) * dt

    tau_c = correlation_time(acf, dt)

    # Error estimation via block averaging
    tau_mean, tau_err, tau_blocks = block_correlation_time(
        interaction_energy, time, n_blocks=5
    )

    # ===== SAVE CSV FILES =====

    np.savetxt(
        "interaction_energy.csv",
        np.column_stack((time, interaction_energy)),
        delimiter=",",
        header="time,interaction_energy",
        comments=""
    )

    np.savetxt(
        "interaction_energy_autocorrelation.csv",
        np.column_stack((lags, acf)),
        delimiter=",",
        header="lag_time,autocorrelation",
        comments=""
    )

    # ===== PLOTS =====

    plt.figure()
    plt.plot(time, interaction_energy)
    plt.xlabel("Time")
    plt.ylabel("Interaction Energy")
    plt.title("Interaction Energy vs Time")
    plt.tight_layout()
    plt.savefig("interaction_energy.png", dpi=300)
    plt.close()

    plt.figure()
    plt.plot(lags, acf)
    plt.axhline(0, linestyle="--", linewidth=0.8)
    plt.xlabel("Lag time")
    plt.ylabel("Normalized autocorrelation")
    plt.title("Autocorrelation of Interaction Energy")
    plt.tight_layout()
    plt.savefig("interaction_energy_autocorrelation.png", dpi=300)
    plt.close()

    # ===== OUTPUT DAT =====

    with open("interaction_autocorrelation_results.dat", "w") as f:
        f.write("Autocorrelation Analysis of Interaction Energy\n\n")

        f.write("Definitions:\n")
        f.write("  E_int(t) = E_Coul(t) + E_LJ(t)\n")
        f.write("  C(tau) = <deltaE(t) deltaE(t+tau)> / <deltaE^2>\n")
        f.write("  tau_c = integral of C(tau) until first zero crossing\n\n")

        f.write(f"Number of frames        = {len(interaction_energy)}\n")
        f.write(f"Time step (dt)          = {dt:.6f}\n")
        f.write(f"Number of blocks        = 5\n\n")

        f.write(f"Correlation time (tau_c)        = {tau_mean:.6f}\n")
        f.write(f"Statistical error (sigma_tau_c) = {tau_err:.6f}\n\n")

        f.write("Block correlation times:\n")
        for i, val in enumerate(tau_blocks, 1):
            f.write(f"  Block {i}: {val:.6f}\n")

    print("Autocorrelation analysis completed successfully.")
    print("Generated files:")
    print("  interaction_energy.png")
    print("  interaction_energy_autocorrelation.png")
    print("  interaction_energy.csv")
    print("  interaction_energy_autocorrelation.csv")
    print("  interaction_autocorrelation_results.dat")


if __name__ == "__main__":
    main()

