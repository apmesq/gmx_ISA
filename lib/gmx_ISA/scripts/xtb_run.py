#usar: python3 xtb_run.py --ligand_charge +1 (NAO REMOVER)

#!/usr/bin/env python3

import os
import sys
import glob
import subprocess
import argparse
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.stats import linregress
from collections import defaultdict

# --- FUNÇÕES AUXILIARES DE PLOTAGEM ---

def plot_smooth_line(x, y, ax, label, color, marker='o'):
    """Plota pontos originais e uma curva suavizada."""
    # Ordenar dados para o spline funcionar
    sorted_indices = np.argsort(x)
    x_sorted = np.array(x)[sorted_indices]
    y_sorted = np.array(y)[sorted_indices]
    
    # Plotar pontos originais
    ax.scatter(x_sorted, y_sorted, color=color, alpha=0.6, label=f"{label} (Points)")
    
    # Criar curva suavizada (Spline)
    if len(x_sorted) > 3:
        try:
            X_Y_Spline = make_interp_spline(x_sorted, y_sorted)
            X_ = np.linspace(x_sorted.min(), x_sorted.max(), 500)
            Y_ = X_Y_Spline(X_)
            ax.plot(X_, Y_, color=color, linewidth=2, label=f"{label} (Smooth)")
        except:
            # Fallback se o spline falhar
            ax.plot(x_sorted, y_sorted, color=color, linewidth=2, label=f"{label}")
    else:
        ax.plot(x_sorted, y_sorted, color=color, linewidth=2, label=f"{label}")

def create_plots(data_rows, output_dir="."):
    """Gera os graficos solicitados."""
    
    # Extrair dados validos (filtrar None)
    frames = []
    e_qm = []
    e_mm = []
    
    for row in data_rows:
        if row['E_Interaction_kJ'] is not None and row['E_Int_MM_kJ'] is not None:
            frames.append(row['Time'])
            e_qm.append(row['E_Interaction_kJ'])
            e_mm.append(row['E_Int_MM_kJ'])
    
    if not frames:
        print("⚠️ Not enough data to generate plots.")
        return

    frames = np.array(frames)
    e_qm = np.array(e_qm)
    e_mm = np.array(e_mm)

    # Configuração geral
    plt.style.use('default') 
    
    # 1. Frame vs E_int QM
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_smooth_line(frames, e_qm, ax, "QM Interaction", "blue")
    ax.set_xlabel("Time (Frame)")
    ax.set_ylabel("Energy (kJ/mol)")
    ax.set_title("QM Interaction Energy vs Time")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.savefig(os.path.join(output_dir, "plot_QM_interaction.png"), dpi=300)
    plt.close(fig)

    # 2. Frame vs E_int MM
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_smooth_line(frames, e_mm, ax, "MM Interaction", "green")
    ax.set_xlabel("Time (Frame)")
    ax.set_ylabel("Energy (kJ/mol)")
    ax.set_title("MM Interaction Energy vs Time")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.savefig(os.path.join(output_dir, "plot_MM_interaction.png"), dpi=300)
    plt.close(fig)

    # 3. Frame vs QM e MM (Juntos)
    fig, ax = plt.subplots(figsize=(8, 6))
    plot_smooth_line(frames, e_qm, ax, "QM", "blue")
    plot_smooth_line(frames, e_mm, ax, "MM", "green")
    ax.set_xlabel("Time (Frame)")
    ax.set_ylabel("Energy (kJ/mol)")
    ax.set_title("QM vs MM Interaction Energy over Time")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.savefig(os.path.join(output_dir, "plot_QM_MM_overlay.png"), dpi=300)
    plt.close(fig)

    # 4. Regressão Linear (QM vs MM)
    fig, ax = plt.subplots(figsize=(8, 6))
    slope, intercept, r_value, p_value, std_err = linregress(e_mm, e_qm)
    r_squared = r_value**2
    
    # Pontos
    ax.scatter(e_mm, e_qm, color='purple', alpha=0.7, label='Data Points')
    
    # Linha de regressão
    fit_line = slope * e_mm + intercept
    ax.plot(e_mm, fit_line, color='black', linestyle='--', linewidth=2, label='Linear Fit')
    
    # Texto com estatísticas
    text_str = f'$R^2 = {r_squared:.3f}$\nSlope = {slope:.3f}\nStd Err = {std_err:.3f}'
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, text_str, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    ax.set_xlabel("MM Interaction Energy (kJ/mol)")
    ax.set_ylabel("QM Interaction Energy (kJ/mol)")
    ax.set_title("Correlation: QM vs MM Energies")
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    fig.savefig(os.path.join(output_dir, "plot_regression.png"), dpi=300)
    plt.close(fig)
    
    print("\n📊 Plots generated successfully.")

# --- FUNÇÕES DO SCRIPT ORIGINAL E NOVAS FUNÇÕES ---

def parse_xvg_file(xvg_path):
    """
    Lê o arquivo coul_lj_sr.xvg.
    Assume formato padrão GROMACS (gmx energy):
    Col 0: Time
    Col 1: Coulomb-SR (ou termo 1)
    Col 2: LJ-SR (ou termo 2)
    Retorna dicionário: {time (int): (coulomb, lj)}
    """
    mm_data = {}
    if not os.path.exists(xvg_path):
        print(f"⚠️ Warning: '{xvg_path}' not found. MM columns will be empty.")
        return mm_data

    print(f"      Reading MM data from {xvg_path}...")
    try:
        with open(xvg_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(('#', '@')):
                    continue
                
                parts = line.split()
                try:
                    # Assume time is first column. Casting to int for matching with frame folder names
                    t = int(float(parts[0])) 
                    # Assume col 1 is Coul and col 2 is LJ (Standard output if requested in that order)
                    coul = float(parts[1])
                    lj = float(parts[2])
                    mm_data[t] = (coul, lj)
                except (ValueError, IndexError):
                    continue
    except Exception as e:
        print(f"❌ Error reading XVG: {e}")
    
    return mm_data

def pdb_to_xyz(pdb_path, xyz_path):
    """Simple PDB to XYZ converter."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip()
                    if not element:
                        element = ''.join([i for i in atom_name if not i.isdigit()])[0]
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append((element, x, y, z))
                except ValueError:
                    continue
    with open(xyz_path, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"Generated from {os.path.basename(pdb_path)}\n")
        for atom in atoms:
            f.write(f"{atom[0]:<4} {atom[1]:12.6f} {atom[2]:12.6f} {atom[3]:12.6f}\n")

def get_smart_protein_charge(pdb_path):
    """Calculates protein charge by inspecting explicit hydrogens."""
    residues = defaultdict(list)
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"): 
                res_name = line[17:20].strip()
                chain = line[21:22].strip()
                res_seq = line[22:26].strip()
                atom_name = line[12:16].strip()
                unique_id = f"{res_name}_{chain}_{res_seq}"
                residues[unique_id].append(atom_name)
    total_charge = 0
    for res_id, atoms in residues.items():
        res_name = res_id.split('_')[0]
        atom_set = set(atoms)
        if res_name == 'ASP':
            total_charge += 0 if any(h in atom_set for h in ['HD2', 'HD1']) else -1
        elif res_name == 'GLU':
            total_charge += 0 if any(h in atom_set for h in ['HE2', 'HE1']) else -1
        elif res_name == 'LYS':
            hz_count = sum(1 for a in atom_set if a.startswith('HZ'))
            total_charge += 1 if hz_count >= 3 else 0 
        elif res_name == 'ARG':
            total_charge += 1
        elif res_name in ['HIS', 'HIE', 'HID', 'HIP']:
            total_charge += 1 if 'HD1' in atom_set and 'HE2' in atom_set else 0
        elif res_name == 'CYS':
            total_charge += 0 if 'HG' in atom_set else -1
        elif res_name == 'TYR':
            total_charge += 0 if 'HH' in atom_set else -1
    return total_charge

def calculate_total_charge(pdb_path, ligand_charge_val):
    filename = os.path.basename(pdb_path)
    ligand_charge_val = int(ligand_charge_val)
    if "ligand" in filename.lower():
        print(f"      > Type: Ligand | Charge: {ligand_charge_val} (User defined)")
        return ligand_charge_val
    protein_charge = get_smart_protein_charge(pdb_path)
    if "complex" in filename.lower():
        total = protein_charge + ligand_charge_val
        print(f"      > Type: Complex | Protein: {protein_charge} + Ligand: {ligand_charge_val} = Total: {total}")
        return total
    print(f"      > Type: Protein | Charge: {protein_charge} (Calculated from atoms)")
    return protein_charge

def run_xtb_calculation(xyz_file, charge, output_dir):
    """Runs GFN2-xTB with robust convergence settings."""
    original_cwd = os.getcwd()
    base_name = os.path.splitext(os.path.basename(xyz_file))[0]
    xtb_work_dir = os.path.join(output_dir, f"xtb_{base_name}")
    if not os.path.exists(xtb_work_dir):
        os.makedirs(xtb_work_dir)
    abs_xyz = os.path.abspath(xyz_file)
    os.chdir(xtb_work_dir)
    cmd = ["xtb", abs_xyz, "--gfn", "2", "--alpb", "water", "--chrg", str(charge), "--etemp", "1500", "--cycles", "1000"]
    print(f"      Running xTB... (Charge: {charge}, etemp: 1500K)")
    try:
        with open(f"{base_name}.out", "w") as outfile:
            subprocess.run(cmd, stdout=outfile, stderr=subprocess.STDOUT, check=True)
        print(f"      ✅ Done.")
    except subprocess.CalledProcessError:
        print(f"      ❌ FAILED. Check output in {xtb_work_dir}")
    except FileNotFoundError:
        print("      ❌ Error: 'xtb' executable not found in PATH.")
        sys.exit(1)
    finally:
        os.chdir(original_cwd)

def extract_energy_from_out(out_file):
    """Extrai a energia total do arquivo de saida do xTB com lógica de fallback robusta."""
    if not os.path.exists(out_file):
        return None
    with open(out_file, 'r') as f:
        lines = f.readlines()
    
    # 1. Tenta pegar a linha final "TOTAL ENERGY"
    for line in reversed(lines):
        if "TOTAL ENERGY" in line and "|" in line:
            parts = line.split()
            try:
                idx = parts.index("ENERGY")
                return float(parts[idx + 1])
            except (ValueError, IndexError):
                pass
    
    # 2. Fallback: Procura convergencia e sobe para achar dados
    convergence_line_idx = -1
    for i, line in enumerate(lines):
        if "*** convergence criteria satisfied after" in line:
            convergence_line_idx = i
            break
    
    if convergence_line_idx != -1:
        for j in range(convergence_line_idx - 1, -1, -1):
            line_content = lines[j].strip()
            if not line_content: continue
            parts = line_content.split()
            if len(parts) >= 2:
                try:
                    int(parts[0]) # Verifica se é iteração
                    return float(parts[1]) # Energia
                except ValueError:
                    continue
    return None

def main():
    parser = argparse.ArgumentParser(description="Batch run GFN2-xTB with Robust Convergence settings.")
    parser.add_argument("--ligand_charge", type=int, default=0, help="Total integer charge of the ligand")
    parser.add_argument("--base_dir", type=str, default="QM_FRAMES_CAPPED", help="Directory containing frame folders")
    args = parser.parse_args()
    
    base_dir = args.base_dir
    if not os.path.exists(base_dir):
        print(f"❌ Error: Directory {base_dir} not found.")
        sys.exit(1)

    frame_dirs = sorted(glob.glob(os.path.join(base_dir, "frame_*")))
    if not frame_dirs:
        print(f"⚠️ No 'frame_*' directories found.")
        sys.exit(0)

    print(f"--- STARTING ROBUST XTB BATCH PROCESSING ---")
    print(f"Global Ligand Charge Set to: {args.ligand_charge}")
    print(f"Settings: --etemp 1500 (Fermi Smearing) active")

    # --- EXECUÇÃO DOS CÁLCULOS (Mantida) ---
    for frame_dir in frame_dirs:
        print(f"\n📂 Processing {os.path.basename(frame_dir)}...")
        pdbs = glob.glob(os.path.join(frame_dir, "*.pdb"))
        pdbs.sort() 
        for pdb in pdbs:
            fname = os.path.basename(pdb)
            print(f"   • File: {fname}")
            total_charge = calculate_total_charge(pdb, args.ligand_charge)
            xyz_name = fname.replace(".pdb", ".xyz")
            xyz_path = os.path.join(frame_dir, xyz_name)
            pdb_to_xyz(pdb, xyz_path)
            run_xtb_calculation(xyz_path, total_charge, frame_dir)

    print("\n--- CALCULATIONS COMPLETED. STARTING DATA EXTRACTION ---")

    # --- PARTE NOVA: Ler arquivo MM (XVG) ---
    xvg_file = "coul_lj_sr.xvg"
    mm_data_map = parse_xvg_file(xvg_file)

    results_csv = "final_energies.csv"
    data_rows = []

    for frame_dir in frame_dirs:
        frame_name = os.path.basename(frame_dir)
        try:
            # Assume formato frame_NUMERO
            time_val = int(frame_name.split('_')[1])
        except:
            # Fallback se nome não for padrao, tenta pegar numeros
            import re
            nums = re.findall(r'\d+', frame_name)
            time_val = int(nums[0]) if nums else 0
            
        energies = {'complex': None, 'ligand': None, 'residues': None}
        xtb_dirs = glob.glob(os.path.join(frame_dir, "xtb_*"))
        
        for x_dir in xtb_dirs:
            dirname = os.path.basename(x_dir)
            if "complex" in dirname.lower(): m_type = 'complex'
            elif "ligand" in dirname.lower(): m_type = 'ligand'
            else: m_type = 'residues'
            
            base_fname = dirname.replace("xtb_", "")
            out_path = os.path.join(x_dir, f"{base_fname}.out")
            val = extract_energy_from_out(out_path)
            if val is not None:
                energies[m_type] = val
        
        # --- CÁLCULOS ENERGIA ---
        e_complex = energies['complex']
        e_ligand = energies['ligand']
        e_residues = energies['residues']
        
        e_int_hartree = None
        e_int_kj = None
        
        if e_complex is not None and e_ligand is not None and e_residues is not None:
            e_int_hartree = e_complex - e_ligand - e_residues
            e_int_kj = e_int_hartree * 2625.5 # Conversão
            
        # --- DADOS MM ---
        mm_coul = None
        mm_lj = None
        e_int_mm_kj = None
        
        if time_val in mm_data_map:
            mm_coul, mm_lj = mm_data_map[time_val]
            e_int_mm_kj = mm_coul + mm_lj

        data_rows.append({
            'Frame': frame_name,
            'Time': time_val,
            'E_Complex_Ha': e_complex,
            'E_Ligand_Ha': e_ligand,
            'E_Residues_Ha': e_residues,
            'E_Interaction_Ha': e_int_hartree,
            'E_Interaction_kJ': e_int_kj,
            'MM_Coulomb_kJ': mm_coul,
            'MM_LJ_kJ': mm_lj,
            'E_Int_MM_kJ': e_int_mm_kj
        })

    # --- ESCREVER CSV ---
    csv_columns = [
        'Frame', 'Time', 
        'E_Complex_Ha', 'E_Ligand_Ha', 'E_Residues_Ha', 
        'E_Interaction_Ha', 'E_Interaction_kJ',
        'MM_Coulomb_kJ', 'MM_LJ_kJ', 'E_Int_MM_kJ'
    ]
    
    try:
        with open(results_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
            writer.writeheader()
            for data in data_rows:
                writer.writerow(data)
        print(f"\n✅ Results saved to: {os.path.abspath(results_csv)}")
    except IOError as e:
        print(f"\n❌ Error writing CSV: {e}")

    # --- GERAR GRÁFICOS ---
    print("\n--- GENERATING PLOTS ---")
    create_plots(data_rows)

    print("\n--- ALL TASKS COMPLETED ---")

if __name__ == "__main__":
    main()
