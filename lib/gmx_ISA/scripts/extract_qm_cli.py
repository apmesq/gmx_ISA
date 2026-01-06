# para usar: python3 extrair_qm_cli.py rerun_com_grupos.tpr center.xtc SUB 6.0 50
# faz o capping automático baseado nas valências de R-C=O e R-NH, sem alterar protonação all atom do gromacs
# necessário que o usuário faça inspeção visual dos cappings antes de prosseguir

#!/usr/bin/env python3

import MDAnalysis as mda
import sys
import os
import warnings
import numpy as np

# Suppress warnings
warnings.filterwarnings("ignore")

def create_hydrogen_caps(universe, selected_residues):
    """
    Detects severed peptide bonds and creates Hydrogen atoms
    to complete valency, based on original geometry.
    """
    new_coords = []
    new_names = []
    new_resnames = []
    new_resids = []
    
    # Set of selected residue IDs for fast lookup
    selected_ids = set(selected_residues.resids)
    
    # Bond lengths (Angstroms)
    dist_nh = 1.01  # N-H (amide)
    dist_ch = 1.10  # C-H (aldehyde/formyl - C-term cap)

    for res in selected_residues:
        # --- 1. CAP AT NITROGEN (N-terminus of the fragment) ---
        # Check if previous residue exists in universe but is NOT in selection
        prev_res_idx = res.resindex - 1
        if prev_res_idx >= 0:
            prev_res = universe.residues[prev_res_idx]
            # If previous residue was NOT selected, there is a cut here
            if prev_res.resid not in selected_ids:
                try:
                    # Atoms involved in the broken bond
                    atom_n = res.atoms.select_atoms("name N")[0]
                    atom_prev_c = prev_res.atoms.select_atoms("name C")[0]
                    
                    # Vector: From previous C towards current N
                    # We want H in the opposite direction: From N back to where C was
                    vec = atom_prev_c.position - atom_n.position
                    norm = np.linalg.norm(vec)
                    unit_vec = vec / norm
                    
                    # New H position: N position + vector * distance
                    h_pos = atom_n.position + (unit_vec * dist_nh)
                    
                    new_coords.append(h_pos)
                    new_names.append("H_cap")
                    new_resnames.append(res.resname)
                    new_resids.append(res.resid)
                except IndexError:
                    pass # Atom not found (e.g., real N-terminus)

        # --- 2. CAP AT CARBON (C-terminus of the fragment) ---
        # Check if next residue exists in universe but is NOT in selection
        next_res_idx = res.resindex + 1
        if next_res_idx < len(universe.residues):
            next_res = universe.residues[next_res_idx]
            # If next residue was NOT selected, there is a cut here
            if next_res.resid not in selected_ids:
                try:
                    # Atoms involved in the broken bond
                    atom_c = res.atoms.select_atoms("name C")[0]
                    atom_next_n = next_res.atoms.select_atoms("name N")[0]
                    
                    # Vector: From current C towards next N
                    vec = atom_next_n.position - atom_c.position
                    norm = np.linalg.norm(vec)
                    unit_vec = vec / norm
                    
                    # New H position: C position + vector * distance
                    h_pos = atom_c.position + (unit_vec * dist_ch)
                    
                    new_coords.append(h_pos)
                    new_names.append("H_cap")
                    new_resnames.append(res.resname)
                    new_resids.append(res.resid)
                except IndexError:
                    pass

    # If no caps were created, return None
    if not new_coords:
        return None

    n_caps = len(new_coords)
    # Create empty Universe for caps (default 1 segment)
    u_caps = mda.Universe.empty(n_atoms=n_caps, n_residues=n_caps, atom_resindex=range(n_caps), trajectory=True)
    
    u_caps.add_TopologyAttr('name', new_names)
    u_caps.add_TopologyAttr('resname', new_resnames)
    u_caps.add_TopologyAttr('resid', new_resids)
    
    # Assign a single generic SEGID to all caps to avoid dimension errors
    u_caps.add_TopologyAttr('segid', ['CAPS']) 

    u_caps.atoms.positions = new_coords
    
    return u_caps.atoms

def main():
    # 1. ARGUMENT CHECK
    if len(sys.argv) != 6:
        print("------------------------------------------------------------------")
        print("❌ ERROR: Incorrect number of arguments.")
        print("USAGE:")
        print("python3 extract_qm_cli.py <tpr> <xtc> <ligand_resname> <cutoff> <stride>")
        print("------------------------------------------------------------------")
        sys.exit(1)

    tpr_file = sys.argv[1]
    xtc_file = sys.argv[2]
    lig_resname = sys.argv[3]
    try:
        cutoff = float(sys.argv[4])
        stride = int(sys.argv[5])
    except ValueError:
        print("❌ ERROR: Cutoff must be float and Stride must be int.")
        sys.exit(1)

    if not os.path.exists(tpr_file) or not os.path.exists(xtc_file):
        print("❌ ERROR: TPR or XTC file not found.")
        sys.exit(1)

    print(f"--- STARTING PROCESSING (WITH H CAPPING) ---")
    
    try:
        u = mda.Universe(tpr_file, xtc_file)
    except Exception as e:
        print(f"❌ ERROR loading trajectory: {e}")
        sys.exit(1)

    base_output_dir = "QM_FRAMES_CAPPED"
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    count = 0

    # 4. LOOP OVER FRAMES
    for ts in u.trajectory[::stride]:
        time_ps = int(ts.time)
        print(f"-> Processing Frame {ts.frame} (t={time_ps}ps)...")

        ligand = u.select_atoms(f"resname {lig_resname}")
        if len(ligand) == 0:
            print(f"   ⚠️ Ligand {lig_resname} not found in this frame.")
            continue

        # Select nearby residues
        nearby_residues = u.select_atoms(f"protein and (around {cutoff} group lig)", lig=ligand).residues
        protein_atoms = nearby_residues.atoms

        # --- GENERATE CAPS ---
        caps_atoms = create_hydrogen_caps(u, nearby_residues)

        # Assemble complex
        if caps_atoms:
            # Merge: Ligand + Protein + New Hydrogens
            complex_atoms = mda.Merge(ligand, protein_atoms, caps_atoms)
            # Pocket only: Protein + New Hydrogens
            pocket_atoms = mda.Merge(protein_atoms, caps_atoms)
            print(f"   + Added {len(caps_atoms)} capping hydrogens.")
        else:
            complex_atoms = mda.Merge(ligand, protein_atoms)
            pocket_atoms = protein_atoms
            print(f"   . No backbone cuts detected.")

        # --- SAVE FILES ---
        frame_dir = os.path.join(base_output_dir, f"frame_{time_ps}ps")
        if not os.path.exists(frame_dir):
            os.makedirs(frame_dir)

        # Updated filenames to include time_ps
        ligand.write(os.path.join(frame_dir, f"ligand_isolated_{time_ps}ps.pdb"))
        pocket_atoms.atoms.write(os.path.join(frame_dir, f"residues_pocket_capped_{time_ps}ps.pdb"))
        complex_atoms.atoms.write(os.path.join(frame_dir, f"complex_capped_{time_ps}ps.pdb"))

        count += 1

    print("------------------------------------------------------------------")
    print(f"✅ SUCCESS! {count} frames processed and capped.")
    print(f"📂 Results saved in folder: ./{base_output_dir}")
    print("------------------------------------------------------------------")

if __name__ == "__main__":
    main()
