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
    rmsf_avg = []
    sasa_avg = []
    step_labels = []

    # ------------------------------
    # LOOP OVER MD STEPS
    # ------------------------------
    for xtc_file in xtc_files:
        xtc_path = os.path.join(ph_path, xtc_file)
        step_name = os.path.splitext(xtc_file)[0]
        print(f"Processing pH {ph}, file {xtc_file}...")

        # ------------------------------
        # RMSF CALCULATION
        # ------------------------------
        rmsf_xvg = os.path.join(ph_path, f"{step_name}_rmsf.xvg")
        cmd_rmsf = [gmx, "rmsf", "-s", tpr_file, "-f", xtc_path, "-o", rmsf_xvg, "-res"]
        subprocess.run(cmd_rmsf, input='1\n', text=True, check=True)

        with open(rmsf_xvg) as f:
            vals = [float(line.split()[1]) for line in f if not line.startswith(('#','@'))]
        rmsf_avg.append(sum(vals)/len(vals))
        
        # ------------------------------
        # SASA CALCULATION
        # ------------------------------
        sasa_xvg = os.path.join(ph_path, f"{step_name}_sasa.xvg")
        total_sasa_xvg = os.path.join(ph_path, f"{step_name}_total_sasa.xvg")
        cmd_sasa = [gmx, "sasa", "-s", tpr_file, "-f", xtc_path,
                    "-o", total_sasa_xvg, "-or", sasa_xvg]
        subprocess.run(cmd_sasa, input='1\n', text=True, check=True)

        with open(sasa_xvg) as f:
            vals = [float(line.split()[1]) for line in f if not line.startswith(('#','@'))]
        sasa_avg.append(sum(vals)/len(vals))
        
        step_labels.append(step_name)
    
    # ------------------------------
    # SAVE DATA TO TEMP FILES FOR GRACE
    # ------------------------------
    temp_rmsf_file = os.path.join(ph_path, "rmsf_tmp.dat")
    temp_sasa_file = os.path.join(ph_path, "sasa_tmp.dat")

    with open(temp_rmsf_file, "w") as f:
        for i, val in enumerate(rmsf_avg, start=1):
            f.write(f"{i} {val}\n")

    with open(temp_sasa_file, "w") as f:
        for i, val in enumerate(sasa_avg, start=1):
            f.write(f"{i} {val}\n")

    # ------------------------------
    # PLOT USING GRACEBAT
    # ------------------------------
    rmsf_plot_png = os.path.join(ph_path, f"rmsf_vs_step_pH{ph}.png")
    sasa_plot_png = os.path.join(ph_path, f"sasa_vs_step_pH{ph}.png")

    subprocess.run([
        "gracebat", "-nxy", temp_rmsf_file,
        "-hdevice", "PNG", "-hardcopy", "-printfile", rmsf_plot_png
    ], check=True)

    subprocess.run([
        "gracebat", "-nxy", temp_sasa_file,
        "-hdevice", "PNG", "-hardcopy", "-printfile", sasa_plot_png
    ], check=True)

    # Remove temporary files
    os.remove(temp_rmsf_file)
    os.remove(temp_sasa_file)

    print(f"Plots saved for pH {ph}: {rmsf_plot_png}, {sasa_plot_png}")

print("All done!")

