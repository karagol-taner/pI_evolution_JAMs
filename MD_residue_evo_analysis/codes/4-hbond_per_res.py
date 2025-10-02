import os
import subprocess
import pandas as pd
import re

os.environ["GMX_MAXBACKUP"] = "-1"

main_folder = "/content/drive/MyDrive/MD/JAM_pH/JAMC/gromacs"
gmx = "/usr/local/gromacs/bin/gmx"
output_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/residue_hbond.csv"

all_results = []

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def xtc_sort_key(filename):
    match = re.search(r'_(\d+)\.xtc$', filename)
    return int(match.group(1)) if match else 0

# Only numeric pH folders
ph_folders = [f for f in os.listdir(main_folder)
              if os.path.isdir(os.path.join(main_folder, f)) and is_float(f)]
ph_folders_sorted = sorted(ph_folders, key=lambda x: float(x))

for ph_folder in ph_folders_sorted:
    ph_path = os.path.join(main_folder, ph_folder)
    
    tpr_files = [f for f in os.listdir(ph_path) if f.endswith(".tpr")]
    xtc_files = [f for f in os.listdir(ph_path) if f.endswith(".xtc")]
    
    if not tpr_files or not xtc_files:
        print(f"No tpr/xtc in {ph_path}")
        continue
    
    tpr_file = os.path.join(ph_path, tpr_files[0])
    xtc_files_sorted = sorted(xtc_files, key=xtc_sort_key)
    
    for xtc_file in xtc_files_sorted:
        xtc_path = os.path.join(ph_path, xtc_file)
        print(f"Processing H-bonds: {ph_folder} - {xtc_file}")
        
        cmd = [
            gmx, "hbond",
            "-s", tpr_file,
            "-f", xtc_path,
            "-num", "hbond.xvg"  # total H-bonds per frame
        ]
        
        # Automatically select protein group (usually 1)
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate(input=b'1\n1\n')  # donor group and acceptor group
        
        if process.returncode != 0:
            print(f"Error in {xtc_file}: {stderr.decode()}")
            continue
        
        # Parse hbond.xvg
        hbond_counts = []
        with open("hbond.xvg") as f:
            for line in f:
                if line.startswith(('#','@')):
                    continue
                parts = line.split()
                hbond_counts.append(float(parts[1]))
        
        avg_hbonds = sum(hbond_counts)/len(hbond_counts) if hbond_counts else 0
        
        all_results.append({
            "pH": ph_folder,
            "xtc_file": xtc_file,
            "avg_hbonds": avg_hbonds
        })

df = pd.DataFrame(all_results)
df.to_csv(output_csv, index=False)
print(f"Saved H-bond results to {output_csv}")

