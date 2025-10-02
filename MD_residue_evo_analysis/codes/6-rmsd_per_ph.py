import os
import subprocess

# ------------------------------
# SETTINGS
# ------------------------------
os.environ["GMX_MAXBACKUP"] = "-1"  # Disable GROMACS backup files
main_folder = "/content/drive/MyDrive/MD/JAM_pH/JAMC/gromacs"
gmx = "/usr/local/gromacs/bin/gmx"

# Ensure gracebat is installed
!apt-get update -qq
!apt-get install -y grace

# ------------------------------
# Detect numeric pH folders only
# ------------------------------
ph_folders = [f for f in os.listdir(main_folder) if os.path.isdir(os.path.join(main_folder, f))]
ph_folders_sorted = sorted([f for f in ph_folders if f.replace('.', '', 1).isdigit()], key=float)
print("Detected pH folders:", ph_folders_sorted)

# ------------------------------
# LOOP OVER PH FOLDERS
# ------------------------------
for ph in ph_folders_sorted:
    ph_path = os.path.join(main_folder, ph)
    tpr_files = sorted([f for f in os.listdir(ph_path) if f.endswith(".tpr")])
    xtc_files = sorted([f for f in os.listdir(ph_path) if f.endswith(".xtc")],
                       key=lambda x: int(''.join(filter(str.isdigit, x))))
    
    if not tpr_files or not xtc_files:
        print(f"No TPR/XTC files found in {ph_path}, skipping...")
        continue
    
    tpr_file = os.path.join(ph_path, tpr_files[0])
    rmsd_avg = []
    step_labels = []

    # ------------------------------
    # LOOP OVER MD STEPS
    # ------------------------------
    for xtc_file in xtc_files:
        xtc_path = os.path.join(ph_path, xtc_file)
        step_name = os.path.splitext(xtc_file)[0]
        print(f"Processing RMSD for pH {ph}, file {xtc_file}...")

        # ------------------------------
        # RMSD CALCULATION
        # ------------------------------
        rmsd_xvg = os.path.join(ph_path, f"{step_name}_rmsd.xvg")
        cmd_rmsd = [gmx, "rms", "-s", tpr_file, "-f", xtc_path, "-o", rmsd_xvg]
        
        # Automatically select backbone group (1) for RMSD reference and fitting
        subprocess.run(cmd_rmsd, input='1\n1\n', text=True, check=True)

        # Read RMSD xvg and calculate average RMSD
        with open(rmsd_xvg) as f:
            vals = [float(line.split()[1]) for line in f if not line.startswith(('#','@'))]
        rmsd_avg.append(sum(vals)/len(vals))
        step_labels.append(step_name)
    
    # ------------------------------
    # SAVE DATA TO TEMP FILE FOR GRACE
    # ------------------------------
    temp_rmsd_file = os.path.join(ph_path, "rmsd_tmp.dat")
    with open(temp_rmsd_file, "w") as f:
        for i, val in enumerate(rmsd_avg, start=1):
            f.write(f"{i} {val}\n")

    # ------------------------------
    # PLOT USING GRACEBAT
    # ------------------------------
    rmsd_plot_png = os.path.join(ph_path, f"rmsd_vs_step_pH{ph}.png")
    subprocess.run([
        "gracebat", "-nxy", temp_rmsd_file,
        "-hdevice", "PNG", "-hardcopy", "-printfile", rmsd_plot_png
    ], check=True)

    # Remove temporary file
    os.remove(temp_rmsd_file)

    print(f"RMSD plot saved for pH {ph}: {rmsd_plot_png}")

print("All RMSD calculations and plots done!")

