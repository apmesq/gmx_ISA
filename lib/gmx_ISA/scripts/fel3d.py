#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

# -----------------------------
# Input arguments
# -----------------------------
if len(sys.argv) != 5:
    print("Usage: python3 fel_pc1_pc2.py <pc1_pc2.xvg> <Temperature_K> <nbins> <sigma_smooth>")
    sys.exit(1)

xvg_file = sys.argv[1]
T = float(sys.argv[2])
bins = int(sys.argv[3])
sigma_smooth = float(sys.argv[4])

# -----------------------------
# Constants
# -----------------------------
kB = 0.008314462618  # kJ/mol/K
Fmax_empty = 5.0    # kJ/mol assigned to empty bins (GROMACS-like behavior)

# -----------------------------
# Read XVG file (PC1, PC2)
# -----------------------------
pc1 = []
pc2 = []

with open(xvg_file) as f:
    for line in f:
        if line.startswith(("#", "@")):
            continue
        cols = line.split()
        pc1.append(float(cols[0]))
        pc2.append(float(cols[1]))

pc1 = np.array(pc1)
pc2 = np.array(pc2)

# -----------------------------
# 2D Histogram (counts)
# -----------------------------
H, xedges, yedges = np.histogram2d(pc1, pc2, bins=bins)

# Probability
P = H / np.sum(H)

# -----------------------------
# Free Energy Landscape
# -----------------------------
F = np.zeros_like(P)
mask = P > 0

F[mask] = -kB * T * np.log(P[mask])
F -= np.min(F[mask])           # normalize minimum to zero
F[~mask] = Fmax_empty          # high energy for empty bins

# Gaussian smoothing (user-defined)
F = gaussian_filter(F, sigma=sigma_smooth)

# -----------------------------
# Grid centers
# -----------------------------
xcenters = 0.5 * (xedges[:-1] + xedges[1:])
ycenters = 0.5 * (yedges[:-1] + yedges[1:])
X, Y = np.meshgrid(xcenters, ycenters, indexing="ij")

# -----------------------------
# Expand energy domain for plotting (±10%)
# -----------------------------
Fmin = np.min(F)
Fmax_plot = np.max(F)
dF = Fmax_plot - Fmin

vmin = Fmin - 0.1 * dF
vmax = Fmax_plot + 0.1 * dF

# -----------------------------
# Save CSV with PC1, PC2, Energy
# -----------------------------
df_points = pd.DataFrame({
    "PC1": X.flatten(),
    "PC2": Y.flatten(),
    "FreeEnergy_kJmol": F.flatten()
})
df_points.to_csv("FEL_PC1_PC2_points.csv", index=False)

df_map = pd.DataFrame(F, index=xcenters, columns=ycenters)
df_map.to_csv("FEL_PC1_PC2_map.csv")

# -----------------------------
# 2D Free Energy Map
# -----------------------------
plt.figure(figsize=(8, 6))
cf = plt.contourf(X, Y, F, levels=50, cmap="viridis",
                  vmin=vmin, vmax=vmax)
plt.colorbar(cf, label="Free Energy (kJ/mol)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("Free Energy Landscape (2D)")
plt.tight_layout()
plt.savefig("FEL_2D.png", dpi=300)
plt.close()

# -----------------------------
# 3D Free Energy Surface
# -----------------------------
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")

surf = ax.plot_surface(
    X, Y, F,
    cmap="viridis",
    linewidth=0,
    antialiased=True
)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("Free Energy (kJ/mol)")
ax.set_zlim(vmin, vmax)
ax.set_title("Free Energy Landscape (3D)")

fig.colorbar(surf, shrink=0.6, label="Free Energy (kJ/mol)")
plt.tight_layout()

# Show before saving
plt.show()

# Save after visualization
fig.savefig("FEL_3D.png", dpi=300)
plt.close(fig)

# -----------------------------
# PC1 vs Energy projection
# -----------------------------
F_pc1 = np.mean(F, axis=1)

plt.figure(figsize=(7, 5))
plt.plot(xcenters, F_pc1)
plt.xlabel("PC1")
plt.ylabel("Free Energy (kJ/mol)")
plt.title("PC1 Free Energy Projection")
plt.tight_layout()
plt.savefig("FEL_PC1_projection.png", dpi=300)
plt.close()

pd.DataFrame({
    "PC1": xcenters,
    "FreeEnergy_kJmol": F_pc1
}).to_csv("FEL_PC1_projection.csv", index=False)

# -----------------------------
# PC2 vs Energy projection
# -----------------------------
F_pc2 = np.mean(F, axis=0)

plt.figure(figsize=(7, 5))
plt.plot(ycenters, F_pc2)
plt.xlabel("PC2")
plt.ylabel("Free Energy (kJ/mol)")
plt.title("PC2 Free Energy Projection")
plt.tight_layout()
plt.savefig("FEL_PC2_projection.png", dpi=300)
plt.close()

pd.DataFrame({
    "PC2": ycenters,
    "FreeEnergy_kJmol": F_pc2
}).to_csv("FEL_PC2_projection.csv", index=False)

print("FEL analysis completed successfully.")
print(f"Parameters used: T = {T} K | bins = {bins} | sigma_smooth = {sigma_smooth}")

