import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# -------------------------------
# File paths
# -------------------------------
rmsf_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/rmsf_per_residue.csv"
sasa_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/residue_sasa.csv"
hbonds_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/residue_hbond.csv"
rmsd_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/rmsd_per_trajectory.csv"

# Output folder for plots
plot_folder = "/content/drive/MyDrive/MD/JAM_pH/JAMC/plots/"
os.makedirs(plot_folder, exist_ok=True)

# -------------------------------
# Load data
# -------------------------------
rmsf_df = pd.read_csv(rmsf_csv)
sasa_df = pd.read_csv(sasa_csv)
hbonds_df = pd.read_csv(hbonds_csv)
rmsd_df = pd.read_csv(rmsd_csv)

# Ensure numeric
for df in [rmsf_df, sasa_df, hbonds_df, rmsd_df]:
    df['pH'] = df['pH'].astype(float)

# -------------------------------
# Compute average per protein per pH
# -------------------------------
avg_rmsf = rmsf_df.groupby('pH')['rmsf'].mean().reset_index()
avg_sasa = sasa_df.groupby('pH')['sasa'].mean().reset_index()
avg_hbonds = hbonds_df.groupby('pH')['avg_hbonds'].mean().reset_index()
avg_rmsd = rmsd_df.groupby('pH')['rmsd_nm'].mean().reset_index()  # average RMSD over time & trajectories

# -------------------------------
# Plot all together in 4 subplots
# -------------------------------
fig, axs = plt.subplots(4, 1, figsize=(8, 16), sharex=True)

sns.lineplot(data=avg_rmsf, x='pH', y='rmsf', marker='o', ax=axs[0])
axs[0].set_title("Average RMSF per Protein vs pH")
axs[0].set_ylabel("RMSF (nm)")
axs[0].grid(True)

sns.lineplot(data=avg_sasa, x='pH', y='sasa', marker='o', color='orange', ax=axs[1])
axs[1].set_title("Average SASA per Protein vs pH")
axs[1].set_ylabel("SASA (nm²)")
axs[1].grid(True)

sns.lineplot(data=avg_hbonds, x='pH', y='avg_hbonds', marker='o', color='green', ax=axs[2])
axs[2].set_title("Average Number of H-bonds per Protein vs pH")
axs[2].set_ylabel("H-bonds")
axs[2].grid(True)

sns.lineplot(data=avg_rmsd, x='pH', y='rmsd_nm', marker='o', color='red', ax=axs[3])
axs[3].set_title("Average RMSD per Protein vs pH")
axs[3].set_xlabel("pH")
axs[3].set_ylabel("RMSD (nm)")
axs[3].grid(True)

plt.tight_layout()
combined_plot_file = os.path.join(plot_folder, "combined_protein_metrics_vs_pH.png")
plt.savefig(combined_plot_file, dpi=300)
plt.show()

print("All plots saved in:", plot_folder)

