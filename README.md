# gmx_ISA (GROMACS Interaction Stability & Analysis)

![Badge](https://img.shields.io/badge/Version-v0.1.0-blue)
![Badge](https://img.shields.io/badge/License-MIT-green)
![Badge](https://img.shields.io/badge/Platform-Linux-lightgrey)

# gmx_ISA

**gmx_ISA** is an automated Bash-Python framework for **quantum-region preparation and analysis from Molecular Dynamics (MD) trajectories**, with a strong focus on **robust QM-ready structure extraction** and **statistical post-processing** for the GROMACS community.

The software provides **automatic QM region extraction** directly from GROMACS trajectories, including **valence-based capping of broken peptide bonds** that preserves the original **all-atom protonation state**. The generated structures are immediately compatible with **semi-empirical quantum calculations using GFN2-xTB**, enabling fast and reproducible energy evaluations on MD-derived configurations.

In addition, `gmx_ISA` delivers **publication-ready plotting**, advanced **statistical analyses of interaction energies** (including autocorrelation, convergence, regression, and stability metrics), and flexible visualization of free-energy landscapes. Finally, the framework also automates a set of **standard GROMACS analyses** (RMSD, RMSF, SASA, radius of gyration, hydrogen bonds), providing a unified and reproducible workflow from trajectory to figures.


---

## ⚡ Key Features

### ⚛ Quantum-Mechanical Workflows (xTB)

* **Automatic QM Region Extraction:**  
  Extraction of ligand-centered QM regions directly from GROMACS trajectories, including all residues within a user-defined cutoff and configurable stride.

* **Valence-Based Automatic Capping:**  
  Automatic capping of severed peptide bonds based on atomic valences, preserving the original **all-atom protonation state** without introducing artificial protonation changes.

* **Semi-Empirical QM Calculations:**  
  Seamless integration with **GFN2-xTB** for fast and reproducible quantum-mechanical energy calculations on MD-derived structures.

---

### 📊 Advanced Analyses & Statistical Post-processing (Python-based)

* **MM Interaction Energy Analysis:**  
  Calculation of Coulomb-SR, Lennard-Jones-SR, and total MM interaction energies, including descriptive statistics.

* **Interaction Stability Analysis (ISA):**  
  Quantification of interaction stability as the ratio between mean interaction energy and its standard deviation.

* **Energy Autocorrelation Analysis:**  
  Autocorrelation functions and integrated autocorrelation times for interaction energies.

* **Free Energy Landscape Visualization:**  
  Generation of 2D and 3D PCA-based free energy landscapes with publication-ready plots.

* **Linear Regression Analysis:**  
  Regression between Coulomb-SR and Lennard-Jones-SR energy components.

* **Diffusion Analysis:**  
  Automatic detection of the linear regime in MSD curves and calculation of diffusion coefficients using the Einstein relation.

* **Free Energy Convergence Analysis:**  
  Cumulative moving average analysis optimized for **gmx_MMPBSA** energy outputs.

---

### ⚙️ Automated GROMACS Analyses

* **Trajectory Pre-processing:**  
  Automatic centering, fitting, and periodic boundary condition (PBC) removal.

* **Structural Analyses:**  
  RMSD (protein backbone, ligand, multiple ligands), RMSF (protein and ligand), radius of gyration, and hydrogen bond analysis.

* **Energetic Analyses:**  
  Total energy, Coulomb short-range, and Lennard-Jones short-range energy terms.

* **Solvation & Spatial Analyses:**  
  Solvent Accessible Surface Area (SASA) and radial distribution functions (RDF) between ligand and solvent.

* **Free Energy Landscapes (GROMACS-based):**  
  PCA-based FEL calculations for protein or ligand systems over user-defined time intervals.


---

## 📦 Dependencies & Installation

### 1. Dependencies

gmx_ISA requires **Python (v3.11)** and external scientific software to operate. The required Python libraries are:

- NumPy  
- Pandas  
- Matplotlib  
- SciPy  
- Scikit-learn  
- MDAnalysis  

In addition, **GROMACS v2023 or newer** is recommended for MD trajectory analysis, and **xTB (v6.7.1)** is required for quantum-mechanical workflows.

---

### 2. Installation

The software can be used locally by cloning the GitHub repository. In this case, the user is responsible for ensuring that all dependencies are properly installed.

```bash
# Clone repository and configure PATH
git clone https://github.com/apmesq/gmx_ISA.git

cd gmx_ISA

chmod +x bin/gmx_ISA

export PATH=$PWD/bin:$PATH

# Create a Conda environment with all dependencies
conda create -n gmx_ISA -c conda-forge python=3.11 numpy pandas matplotlib scipy scikit-learn mdanalysis xtb gromacs=2023
```
---
