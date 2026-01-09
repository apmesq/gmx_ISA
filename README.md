# gmx_ISA (GROMACS Interaction Stability & Analysis)


![Version](https://img.shields.io/badge/Version-v0.5.0-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Plataform](https://img.shields.io/badge/Platform-Linux-lightgrey)
[![Install Test](https://img.shields.io/github/actions/workflow/status/apmesq/gmx_ISA/conda_install_test.yml?label=Install%20Check&logo=github)](https://github.com/apmesq/gmx_ISA/actions/workflows/conda_install_test.yml)
[![Conda Version](https://img.shields.io/conda/vn/apmesq/gmx_isa.svg?color=success&label=Conda%20Version)](https://anaconda.org/apmesq/gmx_isa)
[![Conda Platform](https://img.shields.io/conda/pn/apmesq/gmx_isa.svg?color=blue&label=Conda%20Platform)](https://anaconda.org/apmesq/gmx_isa)
[![Conda Downloads](https://img.shields.io/conda/dn/apmesq/gmx_isa.svg?label=Conda%20Downloads)](https://anaconda.org/apmesq/gmx_isa)


## Table of Contents

- [⚡ Key Features](#-key-features)
  - [⚛ Quantum-Mechanical Workflows (xTB)](#-quantum-mechanical-workflows-xtb)
  - [📊 Advanced Analyses & Statistical Post-processing](#-advanced-analyses--statistical-post-processing-python-based)
  - [⚙️ Automated GROMACS Analyses](#️-automated-gromacs-analyses)
- [📦 Dependencies & Installation](#-dependencies--installation)
  - [1. Dependencies](#1-dependencies)
  - [2. Installation](#2-installation)
- [⏩ Quick Start](#-quick-start)
  - [Interactive configuration using gen_conf.sh](#interactive-configuration-using-gen_confsh)
- [📖 Command Line Helper](#-command-line-helper)
- [📈 Advanced Workflows](#-advanced-workflows)
  - [Basic GROMACS Analyses](#basic-gromacs-analyses)
  - [Additional GROMACS Analyses](#additional-gromacs-analyses-mandatory-flags)
  - [Advanced Analyses Using Python Routines](#advanced-analyses-using-python-routines)
- [QM Workflows with xTB and Automatic Capping](#qm-workflows-with-xtb-and-automatic-capping)
- [✏️ Design Philosophy](#️-design-philosophy)
  - [1. Automation First](#1-automation-first)
  - [2. Bridging Classical and Quantum Worlds](#2-bridging-classical-and-quantum-worlds-qmmm)
  - [3. Modular Architecture](#3-modular-architecture)
  - [4. Operational Safety & Integrity](#4-operational-safety--integrity)
- [💾 Versioning Policy](#-versioning-policy)
- [⚠️ Assumptions & Limitations](#️-assumptions--limitations)
- [📚 Citation](#-citation)
- [⚖️ License](#️-license)
- [🤖 Development & AI Usage Disclaimer](#-development--ai-usage-disclaimer)
- [📬 Contact](#-contact)
- [👥 Contributing](#-contributing)
  - [Developer Guide: Adding New Analysis Flags](#developer-guide-adding-new-analysis-flags)

---

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

#### 2.1 Conda Environment (**RECOMMENDED**)

The easiest way to install `gmx_ISA` and all its dependencies in a single step is via **Conda** (Anaconda/Miniconda).

```bash
# Create a conda env named gmx_ISA
conda create -n gmx_ISA

# Activate the environment
conda activate gmx_ISA

# Automatic installation
conda install -c apmesq -c conda-forge gmx_isa
```
#### 2.2 Local Execution (Git clone)
Alternatively (for development purposes), the software can be used locally by cloning the GitHub repository. In this case, the user is responsible for ensuring that all dependencies are properly installed.

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

To do so, copy the files located in the [`test`](https://github.com/apmesq/gmx_ISA/tree/main/test) directory at the root of the project into your working directory. This directory contains the following files:

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

### Interactive configuration using `gen_conf.sh`

If the user prefers to populate `variables.conf` interactively through the terminal, the auxiliary script `gen_conf.sh` can be invoked. In this case, input files must be provided via command-line flags, and the `-solo` flag should not be used:

```bash
gmx_ISA -s test.tpr -f test.xtc -n test.ndx -en test.edr
```

The following prompt will appear:

```text
Hello! I'm gmx_ISA and I'm glad to help you analyze your molecular dynamics trajectories in few steps, easy and quick! Let's start? (-:
Would you like to use gen_conf interactive bash interface? (yes/no)
```

 By typing `yes`, the index and energy files will be automatically opened, and the user will be asked to manually input the numerical indices corresponding to each system and energy component. For example:

```bash
What is the number of the term 'Ligand2' in your index (*.ndx) file?
> 13
What is the number of the term 'System' in your index (*.ndx) file?
> 0
What is the number for the term 'Coulomb-SR' (it will be followed by '0')?
> 20
What is the number for the term 'Total Energy' (it will be followed by '0')?
> 20
What is the number of the term 'Protein' in your index (*.ndx) file?
> 1
What is the number of the term 'Backbone' in your index (*.ndx) file?
> 4
What is the number of the term 'Ligand' in your index (*.ndx) file?
> 13
What is the number for the term 'Lennard-Jones (LJ-SR)' (it will be followed by '0')?
> 21
```
Once completed, the `variables.conf` file will be written automatically, and the analysis will be executed immediately, generating the same `gmx_ISA-ANALYSIS` output directory as before.

After `variables.conf` has been created, users may run any analysis directly using the `-solo` flag, or simply answer `"no"` when prompted about using the `gen_conf.sh` interface.

---


## 📖 Command Line Helper

You can access the built-in documentation and flag list at any time directly in your terminal:

```bash
gmx_ISA -h
```

<details>
  <summary>Click to see built-in documentation </summary>
  
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

</details>

---

## 📈 Advanced Workflows

For a complete workflow that explores all features of the software, we start from the [`test`](https://github.com/apmesq/gmx_ISA/tree/main/test) directory, using the files provided there, with `variables.conf` already configured and the `-solo` flag enabled to bypass the interactive `gen_conf.sh` interface.

### Basic GROMACS Analyses

As shown previously, we first run the basic analyses whose outputs will serve as inputs for subsequent steps. With the terminal set to the test directory, run:

```bash 
gmx_ISA -solo 
```

This command will generate the `gmx_ISA-ANALYSIS directory`. Regardless of how many analyses are performed, `gmx_ISA` organizes all outputs within the same directory. For organizational clarity, we will create the directory `defaut_gmx_ISA-ANALYSIS` to isolate the analyses:

```bash 
mv gmx_ISA-ANALYSIS defaut_gmx_ISA-ANALYSIS
```

At the end, the following directory structure will be produced:

```text
defaut_gmx_ISA-ANALYSIS/
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

### Additional GROMACS Analyses (Mandatory Flags)

For radius of gyration, SASA, RDF, and FEL calculations, the corresponding flags must be explicitly provided. The user may enable them all at once or individually, as preferred. For FEL calculations, it is mandatory to specify the system (`-sys <protein or ligand>`) and a time interval in picoseconds (`-b`, `-e`):

```bash 
gmx_ISA -solo -gyr -sasa -rdf -fel -sys protein -b 0 -e 1000
```

This step may take approximately 2 minutes due to the computational complexity of the FEL calculation. Once finished, a new `gmx_ISA-ANALYSIS` directory will be created with the outputs. Rename it to `add_fel_gmx_ISA-ANALYSIS`:


```bash 
mv gmx_ISA-ANALYSIS add_fel_gmx_ISA-ANALYSIS
```

The following files should be generated:

``` text
add_fel_gmx_ISA-ANALYSIS/
├── center.xtc
├── ener.xvg
├── enthalpy.xpm
├── entropy.xpm
├── gmx_ISA_debug.log
├── gyrate.xvg
├── prob.xpm
├── prot_eigenval.xvg
├── prot_eigenvec.trr
├── prot_fel.eps
├── prot_fel.xpm
├── prot_pca.xvg
├── rdf_ligand_water.xvg
└── sasa.xvg
```

As before, `.xvg` files follow the standard GROMACS format and can be natively plotted using XMGrace, as commonly done by the community. For the FEL analysis, the two key outputs are `prot_pca.xvg`, containing the PC1 and PC2 values, and `prot_fel.eps`, a 2D image file suitable for quick visual inspection of the results.

### Advanced Analyses Using Python Routines

For complex analyses that are external to GROMACS, each flag must be executed separately.

#### 3D FEL

To visualize the previously generated FEL in 3D and 2D, and to access `.csv` output files for plotting in external software, use the following command pattern:

```bash 
gmx_ISA -fel3d <pc1_pc2.xvg> <Temperature_K> <nbins> <sigma_smooth>
```

Run:

```bash 
gmx_ISA -solo -fel3d add_fel_gmx_ISA-ANALYSIS/prot_pca.xvg 310 50 2
```

Matplotlib will open an interactive graphical window displaying the 3D plot, allowing the user to rotate, zoom, and customize the visualization. After inspection, close the window and the script will terminate automatically, generating a new `gmx_ISA-ANALYSIS` directory, which can be renamed to `fel3d_gmx_ISA-ANALYSIS`:

```bash 
mv gmx_ISA-ANALYSIS fel3d_gmx_ISA-ANALYSIS
```

The following outputs will be generated:

``` text
fel3d_gmx_ISA-ANALYSIS/
├── FEL_2D.png
├── FEL_3D.png
├── FEL_PC1_PC2_map.csv
├── FEL_PC1_PC2_points.csv
├── FEL_PC1_projection.csv
├── FEL_PC1_projection.png
├── FEL_PC2_projection.csv
└── FEL_PC2_projection.png
```

#### Analyses Based on Coulomb_SR and LJ_SR

For MM interaction energy analyses (Coulomb + Lennard-Jones short-range), linear regression of interaction energy, ISA (Interaction Stability Analysis), and interaction energy autocorrelation, the flags `-int`, `-reg`, `-isa`, and `-autocor` are required, respectively. Execute them sequentially: 

```bash 
gmx_ISA -solo -int defaut_gmx_ISA-ANALYSIS/coulomb_SR.xvg defaut_gmx_ISA-ANALYSIS/LJ_SR.xvg
```

```bash 
gmx_ISA -solo -reg defaut_gmx_ISA-ANALYSIS/coulomb_SR.xvg defaut_gmx_ISA-ANALYSIS/LJ_SR.xvg
```

```bash 
gmx_ISA -solo -isa defaut_gmx_ISA-ANALYSIS/coulomb_SR.xvg defaut_gmx_ISA-ANALYSIS/LJ_SR.xvg
```

```bash 
gmx_ISA -solo -autocor defaut_gmx_ISA-ANALYSIS/coulomb_SR.xvg defaut_gmx_ISA-ANALYSIS/LJ_SR.xvg
```

After execution, the `gmx_ISA-ANALYSIS` directory will contain all outputs generated by these flags. For organization, rename it to `cou_lj_based_gmx_ISA-ANALYSIS`, containing the following files:

``` text
cou_lj_based_gmx_ISA-ANALYSIS/
├── autocor_python.log
├── gmx_ISA_debug.log
├── interaction_autocorrelation_results.dat
├── interaction_energy_analysis.dat
├── interaction_energy_autocorrelation.csv
├── interaction_energy_autocorrelation.png
├── interaction_energy_coulomb_vs_lj.png
├── interaction_energy.csv
├── interaction_energy.png
├── interaction_energy_statistics.dat
├── interaction_python.log
├── isa_python.log
└── ISA_results.dat
```

#### Linear Diffusion and Convergence Analyses

For diffusion analysis and calculation of the diffusion coefficient based on the linear region of the MSD, the MSD output must be previously computed. To identify the most linear region for diffusion coefficient estimation, run:

```bash 
gmx_ISA -solo -diff msd_protein.xvg
```

For convergence analysis, `gmx_ISA` has been fully optimized to operate on `gmx_MMPBSA` outputs, automatically detecting the $\Delta G_{complex}$ column, identifying the frames, and computing the cumulative average.

```bash 
gmx_ISA -solo -conv_pbsa test.csv
```

After renaming the `gmx_ISA-ANALYSIS` directory to `diff_conv_gmx_ISA-ANALYSIS`, the following outputs will be obtained:

``` text
diff_conv_gmx_ISA-ANALYSIS/
├── convergence_stats.csv
├── conv_pbsa_python.log
├── cumulative_average.csv
├── diff_python.log
├── gmx_ISA_debug.log
├── msd_protein_best_linear_region.csv
├── msd_protein_diffusion_results.dat
├── msd_protein_msd_all.png
├── msd_protein_msd_best_region.png
├── msd_protein.xvg
├── plot_convergence_check.png
└── plot_cumulative_avg.png
```

## QM Workflows with xTB and Automatic Capping

Finally, the most advanced analyses performed by `gmx_ISA` involve semi-empirical quantum-mechanical calculations using xTB and preparation of QM regions with automatic valence-based capping. A brief demonstration is provided here, as QM calculations require higher computational resources.

To extract the QM region of the protein–ligand system, use the `-extract_qm` flag following this usage pattern:

``` bash
gmx_ISA -extract_qm <ligand_resname> <cutoff_angstroms> <stride>

gmx_ISA -solo -extract_qm SUB 4 100
```

Within a few seconds, a dedicated QM workflow directory (`QM_analysis`) will be created, containing subdirectories corresponding to frames extracted within 4 Å of the ligand (resname SUB), sampled every 100 frames (stride).

The directory structure will be:

``` text 
QM_analysis/
├── QM_FRAMES_CAPPED/
|    ├── frame_0ps/
|    |    ├── complex_capped_0ps.pdb
|    |    ├── ligand_isolated_0ps.pdb
|    |    └── residues_pocket_capped_0ps.pdb
|    ├── frame_2000ps/
|    |    ├── complex_capped_2000ps.pdb
|    |    ├── ligand_isolated_2000ps.pdb
|    |    └── residues_pocket_capped_2000ps.pdb
|    ├── frame_4000ps/
|    |    ├── complex_capped_4000ps.pdb
|    |    ├── ligand_isolated_4000ps.pdb
|    |    └── residues_pocket_capped_4000ps.pdb
└── center.xtc
```

To compute the interaction energy between the ligand and directly contacting amino acid residues using xTB, `gmx_ISA` processes each capped `.pdb` file, converts it to .`xyz`, submits each component of each frame to xTB, organizes the resulting energy values, and computes the QM interaction energy per frame. Run:

``` bash 
gmx_ISA -xtb_run --ligand_charge <int>

gmx_ISA -solo -xtb_run --ligand_charge +1
```

The output is a `.csv` file containing xTB energies for the complex, ligand, residues, QM ligand–residue interaction energy, MM Coulomb energy, MM Lennard-Jones energy, and total MM interaction energy for each analyzed frame. The QM interaction energy for frame $i$ is computed as:

$$

E_i^{int,\,QM} = E_i^{complex, \,QM} - (E_i^{residues,\,QM} + E_i^{ligand,\,QM})

$$

The MM interaction energy for frame $i$ is computed as:

$$

E_i^{int,\,MM} = E_i^{Coulomb-SR} + E_i^{LJ-SR}

$$

Additionally, high-quality PNG images are generated for *frame vs. MM energy*, *frame vs. QM energy*, *frame vs. QM/MM comparison*, and *linear regression of QM versus MM interaction energies*. The final directory will contain:

``` text
QM_analysis/
├── center.xtc
├── coul_lj_sr.xvg
├── final_energies.csv
├── plot_MM_interaction.png
├── plot_QM_interaction.png
├── plot_QM_MM_overlay.png
├── plot_regression.png
└── QM_FRAMES_CAPPED/
```

### ***All outputs produced in this advanced workflow are available in the [`output_examples`](https://github.com/apmesq/gmx_ISA-examples-docs/tree/main/outputs_examples) directory in the auxiliary repository***

---


## ✏️ Design Philosophy

**gmx_ISA** was built to bridge the gap between raw molecular dynamics trajectories and actionable scientific insight. The development is guided by four core principles designed to help researchers in computational chemistry.

### 1. Automation First
I believe researchers should spend their time interpreting data, not wrestling with repetitive command-line flags.
* **Goal:** Reduce the time from "simulation finished" to "publication-ready plot".
* **Implementation:** The tool automates GROMACS workflows (like RMSD, RMSF, H-Bonds, and Energy calculations) into single-flag operations, ensuring consistency across different projects.

### 2. Bridging Classical and Quantum Worlds (QM/MM)
Modern computational chemistry often requires validating force-field approximations with quantum mechanical data. **gmx_ISA** democratizes access to hybrid workflows.
* **xTB Integration:** Priority in the integration of semi-empirical methods (GFN2-xTB) directly into the end-state MD pipeline. This allows users to seamlessly transition from classical trajectories to quantum property calculations.
* **Smart Extraction & Capping:** Preparing a system for QM should not be a complete manual burden. The software features intelligent algorithms to extract specific regions (QM clusters) and **automatically cap broken peptide bonds**. This ensures that truncated protein segments remain chemically valid and stable during SCF convergence, without requiring manual editing of PDB files.

### 3. Modular Architecture
The software acts as a robust "hub" rather than a monolithic block.
* **Extensibility:** `gmx_ISA` is designed as a wrapper that orchestrates native GROMACS tools and external Python modules.
* **Plug-and-Play:** New analysis methods can be added as isolated modules (flags) without destabilizing the core kernel. This separation allows for the continuous integration of advanced statistical tools (like PCA, clustering, or regressions).

### 4. Operational Safety & Integrity
Speed should never compromise data integrity.
* **Exclusive Execution Locks:** To prevent memory corruption or file conflicts, the system enforces a strict "one-at-a-time" policy for heavy Python-based routines.
* **Pre-Flight Checks:** The software validates the environment, directory structure, and input files before executing any logic, preventing "silent failures" midway through batch processing.

---

##  💾 Versioning Policy

`gmx_ISA` is the result of approximately one year of independent development. The current version available in this repository (**v0.5.0**) has been manually tested by the developer and is considered internally stable, having passed extensive internal validation.

The **v1.0.0** release will be issued after external validation by users willing to contribute to the project and upon publication of the software in a peer-reviewed scientific journal.

Based on these principles, `gmx_ISA` follows **Semantic Versioning (SemVer)**, adopting the following versioning policy:

- **v0.x**: Development phase. Features and internal design may change; the software is internally stable and open to external user validation and contributions.
- **v1.0.0**: First stable public release. The API and analysis workflows are considered stable and suitable for production and publication.
- **v1.x**: Stable releases with minor updates, performance improvements, and new features that preserve backward compatibility.
- **v2.0.0**: Major release introducing breaking changes to the API, workflows, or data formats.

All releases include a changelog detailing new features, bug fixes, and breaking changes.

---

## ⚠️ Assumptions & Limitations

`gmx_ISA` provides fast and reproducible post-processing workflows for MD trajectories. Users should be aware of the following methodological assumptions and limitations.

### Semi-empirical Quantum Calculations

All QM calculations rely on the semi-empirical GFN2-xTB method. While computationally efficient and well-suited for high-throughput analyses, **absolute energies are less accurate than DFT or ab initio approaches**. Results should therefore be interpreted **comparatively rather than quantitatively**.

### Post-simulation QM Analysis

QM calculations are performed **a posteriori** on MD snapshots. This approach is suitable for interaction analysis and force-field validation, but **does not replace fully coupled, on-the-fly QM/MM simulations**, which account for dynamic polarization effects.

### QM Region Definition and Capping

QM regions are generated automatically using valence-based hydrogen capping. Users must **inspect all generated fragments**, verifying hydrogen placement and chemical consistency.  For small QM regions or high-accuracy studies, **ACE/NME-like capping strategies may be preferable**.

### Scope of the QM Workflow

The current QM workflow assumes a **single ligand–protein system** per analysis. More complex systems (multiple ligands, metal centers, supramolecular assemblies) are outside the present scope but may be supported in future versions.

### Dependence on MD Simulation Quality

All results depend on the **quality of the underlying MD simulation**, including equilibration, sampling, force-field choice, and ligand parameterization. Post-processing cannot compensate for poorly converged or inadequately parameterized simulations.

### Statistical Assumptions

Statistical analyses assume:
- Stationary trajectories over the analyzed window
- Adequate sampling
- Correct handling of correlated time series

Violations of these assumptions may lead to misleading results.

### Diffusion Coefficient Estimation

Diffusion coefficients are computed from the **mean squared displacement (MSD)** using the Einstein relation, assuming normal diffusion:

$$
D = \lim_{t \to \infty} \frac{1}{2N} \frac{d}{dt}
\left\langle \left| \mathbf{r}(t) - \mathbf{r}(0) \right|^2 \right\rangle
\approx \frac{1}{6}\frac{d}{dt}MSD(t)
$$

No finite-size corrections are applied. In particular, the **Yeh–Hummer correction** is not included. For quantitative comparison with experimental data, the following correction should be applied:

$$
D_{\infty} \approx \frac{1}{6}\frac{d}{dt}MSD(t) + \frac{k_B T \, \xi}{6 \pi \eta L}
$$

As a result:
- Diffusion coefficients may be systematically underestimated without proper correction
- Values should be interpreted **comparatively** with multiple simulation samples
- Users must ensure a well-defined linear MSD regime.

---


## 📚 Citation
If you use **gmx_ISA** in your research, please cite the software (DOI below) and the underlying engines used for your calculations:

**1. gmx_ISA Software:**
> Mesquita, A. P. L. (2026). gmx_ISA: Automated Framework for MD-QM/MM Analysis.

**2. If you used MD Analysis features, please also cite GROMACS:**
> Abraham, M. J., et al. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*, 1-2, 19-25.

**3. If you used Quantum Workflows (`-xtb_run`), please also cite xTB:**
> Bannwarth, C., et al. (2019). GFN2-xTB—An accurate and broadly parametrized self-consistent tight-binding quantum chemical method. *J. Chem. Theory Comput.*, 15(3), 1652-1671.

---

## ⚖️ License

**gmx_ISA** is open-source software licensed under the **MIT License**.

You are free to use, modify, and distribute this software, provided that the original copyright notice and permission notice are included.

For full legal details, please refer to the [LICENSE](https://github.com/apmesq/gmx_ISA/blob/main/LICENSE) file located in the root directory of this repository.

---

## 🤖 Development & AI Usage Disclaimer

**gmx_ISA** is an independent, authorial project conceived, developed, and validated by **Antônio P. L. de Mesquita**, without external financial support or a dedicated development team.

To accelerate development and ensure code robustness, Large Language Models (LLMs) were utilized strictly as **coding assistants**. Their role was limited to:
- Reviewing and optimizing Bash and Python scripts (especially syntax and error handling).
- Assisting in debugging and resolving compilation errors.
- Improving code readability and documentation clarity.

**Crucially, algorithms, scientific logic, statistical methodologies, and validation tests are original intellectual property of the author and were manually verified.** The AI tools served to refine the implementation, not to generate the scientific content.

---

## 📬 Contact

**Main Developer:** Antônio P. L. de Mesquita (PhD student in computational chemistry at Universidade Federal de Lavras - Brazil)   
**e-mail:** antoniolemos.mesquita@gmail.com

---

## 👥 Contributing

Contributions are welcome! If you find a bug or have a feature request, please [open an issue](https://github.com/apmesq/gmx_ISA/issues) on GitHub.

### Developer Guide: Adding New Analysis Flags

`gmx_ISA` mother script is a monolitic bash script. To add a new analysis feature (whether a Bash one-liner or a complex Python module), follow this 4-step workflow to ensure integration with the CLI and safety checks.

#### 1. Register the Variable

Initialize your control variable in Section 5 (Variable Parsing), alongside the other default boolean states.

``` bash
# Inside variables initialization
DO_MY_FEATURE=false
```

#### 2. Add Flag Parsing

Add your flag to the while loop in Section 5. You must set your control variable to true, mark `CUSTOM_FLAG_USED=true`, and handle argument shifting (shift) if your flag takes parameters.

``` bash
-my_feature)
        DO_MY_FEATURE=true
        CUSTOM_FLAG_USED=true
        # If your flag takes an argument (e.g., input file):
        # my_input="$2"
        # shift 2
        shift 
        ;;
```

#### 3. Safety Check (Python Modules Only)

**Crucial:** If your new flag triggers an external Python script, you must add it to the Exclusive Execution Lock block to prevent memory conflicts.

``` bash
# Inside the py_analysis_count block
if [ "$DO_MY_FEATURE" = true ]; then ((py_analysis_count++)); fi
```

#### 4. Implement Execution Logic

Add your execution block in Section 7. Use the built-in `progress_bar_real_time` function to maintain the UI standard.

**For Bash commands:**

``` bash 
if [ "$DO_MY_FEATURE" = true ]; then
    echo "Running my new analysis..."
    progress_bar_real_time "gmx analyze -f ... -o output.xvg"
fi
```

**For Python scripts:** Ensure your script is located in [`$PYTHON_DIR`](https://github.com/apmesq/gmx_ISA/tree/main/lib/gmx_ISA/scripts) and handle input validation before execution.

``` bash
if [ "$DO_MY_FEATURE" = true ]; then
    # Validate inputs first
    if [ ! -f "$my_input" ]; then echo "Error: File not found"; exit 1; fi
    
    # Execute with progress bar
    SCRIPT_PATH="$PYTHON_DIR/my_script.py"
    progress_bar_real_time "python3 \"$SCRIPT_PATH\" \"$my_input\" > my_log.log"
fi
```

#### 5. File Cleanup

Finally, ensure all generated files (.xvg, .dat, .log) are moved to the analysis folder in Section 8, with simple `mv` commands.

---

