#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import io
import sys

# ===========================================================
# CONFIGURATION: INPUT HANDLING
# ===========================================================
# Check if input file was provided via command line
if len(sys.argv) < 2:
    print("Error: No input file provided.")
    print("Usage: python3 conv_pbsa.py <input_file.csv>")
    sys.exit(1)

INPUT_FILE = sys.argv[1]

# Convergence threshold (kcal/mol)
THRESHOLD = 0.5 

# Output filenames
OUT_CSV_CUMULATIVE = "cumulative_average.csv"
OUT_CSV_INFO = "convergence_stats.csv"
OUT_PLOT_AVG = "plot_cumulative_avg.png"
OUT_PLOT_CHECK = "plot_convergence_check.png"

# ===========================================================
# FUNCTION: Parse gmx_MMPBSA Specific Block
# ===========================================================
def parse_mmpbsa_csv(filename):
    """
    Reads the file line by line to find the 'Delta Energy Terms' block,
    then parses the subsequent CSV data.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    start_index = -1
    for i, line in enumerate(lines):
        if "Delta Energy Terms" in line:
            # The header is usually the line immediately after this title
            start_index = i + 1 
            break
    
    if start_index == -1:
        raise ValueError("Block 'Delta Energy Terms' not found in the file.")

    # Create a string buffer from the relevant lines
    # We join all lines starting from the header line found
    data_str = "".join(lines[start_index:])
    
    # Read into pandas
    try:
        df = pd.read_csv(io.StringIO(data_str), sep=",")
        
        # Clean column names (strip whitespace if any)
        df.columns = df.columns.str.strip()
        
        return df
    except Exception as e:
        raise RuntimeError(f"Error parsing CSV data block: {e}")

# ===========================================================
# 1. Load Data and Calculate Statistics
# ===========================================================
print(f"Loading data from: {INPUT_FILE}")

if not os.path.exists(INPUT_FILE):
    print(f"Error: File '{INPUT_FILE}' not found.")
    sys.exit(1)

try:
    df = parse_mmpbsa_csv(INPUT_FILE)
    
    # Check for the specific columns from gmx_MMPBSA
    # Usually 'Frame #' for X and 'TOTAL' for Delta G
    required_cols = ["Frame #", "TOTAL"]
    
    if not all(col in df.columns for col in required_cols):
        print(f"Error: Columns {required_cols} not found in the parsed block.")
        print(f"Found columns: {list(df.columns)}")
        sys.exit(1)
        
except Exception as e:
    print(f"Processing Error: {e}")
    sys.exit(1)

# Extract Data
G = df["TOTAL"].values
frames = df["Frame #"].values

# Calculate Cumulative Moving Average
# Formula: (Sum of elements up to i) / (count up to i)
cumulative_avg = np.cumsum(G) / np.arange(1, len(G) + 1)

# Calculate Absolute Difference between consecutive averages (Fluctuation)
diff_avg = np.abs(np.diff(cumulative_avg))

# Identify Convergence Point (CAn)
# Find indices where the fluctuation is below the threshold
converged_indices = np.where(diff_avg < THRESHOLD)[0]

# If indices exist, the first one is the convergence point (+1 because diff reduces length by 1)
CAn = converged_indices[0] + 1 if len(converged_indices) > 0 else None

if CAn:
    print(f"Convergence detected at Frame: {CAn}")
else:
    print("Convergence threshold not reached.")

# ===========================================================
# 2. Save Numerical Results to CSV
# ===========================================================
# Save cumulative average data
pd.DataFrame({
    "Frame": frames, 
    "CumulativeAverage": cumulative_avg
}).to_csv(OUT_CSV_CUMULATIVE, index=False)

# Save convergence statistics info
if CAn:
    pd.DataFrame({
        "Threshold": [THRESHOLD], 
        "CAn_Frame": [CAn],
        "Final_Avg_Energy": [cumulative_avg[-1]]
    }).to_csv(OUT_CSV_INFO, index=False)

print("Numerical results saved to CSV files.")

# ===========================================================
# 3. Plotting: Graph 1 - Cumulative Average
# ===========================================================
plt.figure(figsize=(8, 6))

plt.plot(frames, cumulative_avg, color='#1f77b4', linewidth=2, label="Cumulative Avg $\Delta G$")

if CAn:
    plt.axvline(CAn, color="red", linestyle="--", alpha=0.7, label=f"CAn (Frame {CAn})")

plt.xlabel("Frame", fontsize=12)
plt.ylabel("Average $\Delta G$ (kcal/mol)", fontsize=12)
plt.title("Cumulative Convergence Analysis (CAn)", fontsize=14, fontweight='bold')
plt.legend(loc='best')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig(OUT_PLOT_AVG, dpi=300, bbox_inches='tight')
plt.close() 
print(f"Plot 1 saved: {OUT_PLOT_AVG}")

# ===========================================================
# 4. Plotting: Graph 2 - Fluctuation / Difference
# ===========================================================
plt.figure(figsize=(8, 6))

# Note: frames[1:] because diff array is 1 element shorter
plt.plot(frames[1:], diff_avg, color='#ff7f0e', linewidth=1.5, label="|$\Delta$ Moving Avg|")

plt.axhline(THRESHOLD, color="black", linestyle=":", label=f"Threshold ({THRESHOLD})")

if CAn:
    plt.axvline(CAn, color="red", linestyle="--", alpha=0.7, label=f"Convergence Point")

plt.xlabel("Frame", fontsize=12)
plt.ylabel("Change in Average (Fluctuation)", fontsize=12)
plt.title("Convergence Stability Check", fontsize=14, fontweight='bold')
plt.yscale('log') 
plt.legend(loc='best')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig(OUT_PLOT_CHECK, dpi=300, bbox_inches='tight')
plt.close()
print(f"Plot 2 saved: {OUT_PLOT_CHECK}")

print("Analysis complete.")
