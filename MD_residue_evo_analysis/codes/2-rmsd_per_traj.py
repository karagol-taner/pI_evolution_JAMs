import os
import subprocess
import pandas as pd
import re

# -------------------------------
# Settings
# -------------------------------
os.environ["GMX_MAXBACKUP"] = "-1"
main_folder = "/content/drive/MyDrive/MD/JAM_pH/JAMC/gromacs"
gmx = "/usr/local/gromacs/bin/gmx"
output_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/rmsd_per_trajectory.csv"

# -------------------------------
# Helper functions
# -------------------------------
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def xtc_sort_key(filename):
    match = re.search(r'_(\d+)\.xtc$', filename)
    return int(match.group(1)) if match else 0

# -------------------------------
# Find numeric pH folders
# -------------------------------
ph_folders = [f for f in os.listdir(main_folder)
              if os.path.isdir(os.path.join(main_folder, f)) and is_float(f)]
ph_folders_sorted = sorted(ph_folders, key=lambda x: float(x))

all_results = []

# -------------------------------
# Loop over pH folders
# -------------------------------
for ph_folder in ph_folders_sorted:
    ph_path = os.path.join(main_folder, ph_folder)
    
    tpr_files = [f for f in os.listdir(ph_path) if f.endswith(".tpr")]
    xtc_files = [f for f in os.listdir(ph_path) if f.endswith(".xtc")]
    
    if not tpr_files or not xtc_files:
        print(f"No TPR/XTC files in {ph_path}, skipping...")
        continue
    
    tpr_file = os.path.join(ph_path, tpr_files[0])
    xtc_files_sorted = sorted(xtc_files, key=xtc_sort_key)
    
    for xtc_file in xtc_files_sorted:
        xtc_path = os.path.join(ph_path, xtc_file)
        print(f"Processing RMSD: pH {ph_folder}, file {xtc_file}...")
        
        rmsd_xvg = os.path.join(ph_path, f"{xtc_file}_rmsd.xvg")
        
        cmd_rmsd = [
            gmx, "rms",  # RMS calculation
            "-s", tpr_file,
            "-f", xtc_path,
            "-o", rmsd_xvg
        ]
        
        # Automatically select "Backbone" (1) for RMSD
        process = subprocess.Popen(cmd_rmsd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(input=b'1\n1\n')  # Reference + fitting group
        if process.returncode != 0:
            print(f"Error in {xtc_file}: {stderr.decode()}")
            continue
        
        # Parse RMSD xvg
        with open(rmsd_xvg) as f:
            for line in f:
                if line.startswith(('#','@')):
                    continue
                parts = line.split()
                time = float(parts[0])
                rmsd = float(parts[1])
                all_results.append({
                    "pH": ph_folder,
                    "xtc_file": xtc_file,
                    "time_ns": time,
                    "rmsd_nm": rmsd
                })

# -------------------------------
# Save CSV
# -------------------------------
df = pd.DataFrame(all_results)
df.to_csv(output_csv, index=False)
print(f"Saved RMSD results to {output_csv}")

