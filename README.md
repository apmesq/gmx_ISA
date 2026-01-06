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

## ⏩ Quick Start

For a quick execution, after installing `gmx_ISA` and its dependencies, users can directly run the **default analyses performed via GROMACS**. These include: trajectory centering and periodic boundary condition (PBC) removal; calculation of short-range Coulomb and Lennard-Jones energies; total energy (kinetic + potential); protein RMSD; RMSD of up to two ligands present in the trajectory; RMSF of protein and ligand; and the number of hydrogen bonds between protein and ligand.

To do so, copy the files located in the `test` directory at the root of the project into your working directory. This directory contains the following files:

- `test.edr`
- `test.ndx`
- `test.tpr`
- `test.xtc`
- `variables.conf`

Then, execute:

```bash
# Activate the conda environment
conda activate gmx_ISA

# Run the software
gmx_ISA -solo
```

This command will automatically generate the gmx_ISA-ANALYSIS directory, containing the 11 expected output files:

```text

gmx_ISA-ANALYSIS/
├── center.xtc
├── coulomb_SR.xvg
├── gmx_ISA_debug.log
├── LIG-LIG-2_rmsd.xvg
├── LIG-LIG_rmsd.xvg
├── LJ_SR.xvg
├── numLigH.xvg
├── PROT-PROT_rmsd.xvg
├── rmsf_LIG.xvg
├── rmsf_PROT.xvg
└── total_energy.xvg

``` 

The -solo flag forces the analysis to be executed directly. The configuration file variables.conf is already properly set for the test system and reflects the numerical indices defined in the simulation index file (test.ndx). The configuration file follows the structure below:

```text

# Please, modify here according to your *.ndx file only if you won't use gen_conf interactive interface.
Protein="1"
Ligand="13"
System="0"
Backbone="4"
Ligand2="13"
TotalE="2 0"
CoulombSR="20 0"
LJSR="21 0"
# Please, modify the sections bellow according to the output files (*.tpr, *.xtc, *.ndx, *.edr) of your MD simulation.
tpr="./test.tpr"
xtc="./test.xtc"
xtcCENTER="center.xtc" # Please, don't change this line.
ndx="./test.ndx"
edr="./test.edr"

```

## 📖 Command Line Helper

You can access the built-in documentation and flag list at any time directly in your terminal:

```bash
gmx_ISA -h
```

```text
SYNOPSIS

gmx_ISA [-s [<.tpr>]] [-f [<.xtc>]] [-n [<.ndx>]]
        [-en [<.edr>]] [-center] [-rms_lig] 
        [-rms_lig2] [-rms_prot] [-rmsf_prot] 
        [-rmsf_lig][-hbond] [-etot] [-ecoul] 
        [-elj] [-gyr] [-sasa] [-fel] [-sys <prot,lig>] 
        [-e <time (ps)>] [-b <time (ps)>] [-rdf] [-solo]
        [-int [<Coul.xvg>] [<LJ.xvg>]] [-fel3d 
        [<PC1_PC2.xvg>] <temp> <nbins> <smooth>]
        [-reg [<Coul.xvg>] [<LJ.xvg>]] 
        [-isa [<Coul.xvg>] [<LJ.xvg>]] 
        [-autocor [<Coul.xvg>] [<LJ.xvg>]] 
        [-diff [<msd.xvg>]] [-conv_pbsa [<.csv>]] 
        [-extract_qm <lig_resname> <cutoff_angst> <stride>]
        [-run_xtb --ligand_charge <int>]

OPTIONS

Options to specify input files:

    -s                Input <.tpr> file;
    -f                Input <.xtc> file;
    -n                Input <.ndx> file;
    -en               Input <.edr> file.

Default analyses. When no flags are provided, the following
analyses are performed by default:

    -center           Centralize trajectory and remove PBC;
    -rms_lig          Calculates de RMSD of your ligand;
    -rms_lig2         Calculates the RMSD of a second ligand;
    -rms_prot         Calculates the RMSD of the Protein;
    -rmsf_prot        Calculates the RMSD of the Protein;
    -rmsf_lig         Calculates the RMSD of the ligand;
    -hbond            Calculates the number of hydrogen bonds;
    -etot             Calculates the total energy of the system;
    -ecoul            Calculates the Coulomb Short-Range energy;
    -elj              Calculates the LJ Short-Range energy.

Additional analyses performed using GROMACS

    -gyr              Protein's Radius of Gyration;
    -sasa             Protein's Solvent Accessible Surface Area;
    -fel              Free Energy Landscapes through PC1 and PC2;
    -sys              System to perform FEL <ligand, protein>;
    -b                Initial time of the trajectory in picosseconds;
    -e                Final time of the trajectory in picosseconds;
    -rdf              Calculates RDF between ligand and water (OW);
    -solo             Use this flag if you wouldn't like to use gen_conf.
    
External analyses implemented with Python modules. 
Please run one flag at a time.

    -int              Calculates MM interaction energy and statistics;
    -fel3d            Plots 3D and 2D free energy landscape from PC1/2;
    -reg	          Performs linear regression for interaction energy;
    -isa              Calculates the ISA index;
    -autocor          Autocorrelation of the interaction energy;
    -diff             Detects MSD's region with the highest linear;
    -conv_pbsa        Calculates the Cumulative Moving Average.
    
Options for QM calculations using xTB 

    -extract_qm       Extracts the QM region from the <.xtc> file;
    -xtb_run          Semi-empirical quantum calculation with xTB;
    --ligand_charge   Specify the ligand charge.
```

---
