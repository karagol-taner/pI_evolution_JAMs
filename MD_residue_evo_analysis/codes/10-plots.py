import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ============================
# File paths
# ============================
jama_rmsf   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMA/rmsf_per_residue-JAMA.csv"
jama_sasa   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMA/residue_sasa-JAMA.csv"
jama_hbonds = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMA/residue_hbond-JAMA.csv"
jama_rmsd   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMA/rmsd_per_trajectory-JAMA.csv"

jamb_rmsf   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMB/rmsf_per_residue-JAMB.csv"
jamb_sasa   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMB/residue_sasa-JAMB.csv"
jamb_hbonds = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMB/residue_hbond-JAMB.csv"
jamb_rmsd   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMB/rmsd_per_trajectory-JAMB.csv"

jamc_rmsf   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMC/rmsf_per_residue-JAMC.csv"
jamc_sasa   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMC/residue_sasa-JAMC.csv"
jamc_hbonds = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMC/residue_hbond-JAMC.csv"
jamc_rmsd   = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMC/rmsd_per_trajectory-JAMC.csv"

plot_folder = "/Users/tanerkaragol/Desktop/JAM/MD/analysis/plots/"
os.makedirs(plot_folder, exist_ok=True)

# ============================
# Load & label helper
# ============================
def load_and_label(path, protein):
    df = pd.read_csv(path)
    df["pH"] = df["pH"].astype(float)
    df["protein"] = protein
    return df

# ============================
# Load & combine data
# ============================
rmsf_df = pd.concat([load_and_label(jama_rmsf, "JAMA"),
                     load_and_label(jamb_rmsf, "JAMB"),
                     load_and_label(jamc_rmsf, "JAMC")], ignore_index=True)

sasa_df = pd.concat([load_and_label(jama_sasa, "JAMA"),
                     load_and_label(jamb_sasa, "JAMB"),
                     load_and_label(jamc_sasa, "JAMC")], ignore_index=True)

hbonds_df = pd.concat([load_and_label(jama_hbonds, "JAMA"),
                       load_and_label(jamb_hbonds, "JAMB"),
                       load_and_label(jamc_hbonds, "JAMC")], ignore_index=True)

rmsd_df = pd.concat([load_and_label(jama_rmsd, "JAMA"),
                     load_and_label(jamb_rmsd, "JAMB"),
                     load_and_label(jamc_rmsd, "JAMC")], ignore_index=True)

# ============================
# Compute averages
# ============================
avg_rmsf   = rmsf_df.groupby(["protein", "pH"])["rmsf"].mean().reset_index()
avg_sasa   = sasa_df.groupby(["protein", "pH"])["sasa"].mean().reset_index()
avg_hbonds = hbonds_df.groupby(["protein", "pH"])["avg_hbonds"].mean().reset_index()
avg_rmsd   = rmsd_df.groupby(["protein", "pH"])["rmsd_nm"].mean().reset_index()

# ============================
# Custom color palette
# ============================
custom_palette = {"JAMA": "#1b9e77", "JAMB": "#d95f02", "JAMC": "#7570b3"}

# ============================
# Plot 2x2 grid with RMSD first
# ============================
fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

# --- RMSD ---
sns.lineplot(data=avg_rmsd, x="pH", y="rmsd_nm", hue="protein",
             marker="o", palette=custom_palette, ax=axs[0])
axs[0].set_title("Average RMSD vs pH", fontweight='bold')
axs[0].set_ylabel("RMSD (nm)", fontweight='bold')
axs[0].set_xlabel("pH", fontweight='bold')
axs[0].grid(True, linestyle="--")
for tick in axs[0].get_xticklabels() + axs[0].get_yticklabels():
    tick.set_fontweight('bold')
leg = axs[0].legend(title="Protein")
leg.get_title().set_fontweight('bold')

# --- RMSF ---
sns.lineplot(data=avg_rmsf, x="pH", y="rmsf", hue="protein",
             marker="o", palette=custom_palette, ax=axs[1])
axs[1].set_title("Average RMSF vs pH", fontweight='bold')
axs[1].set_ylabel("RMSF (nm)", fontweight='bold')
axs[1].set_xlabel("pH", fontweight='bold')
axs[1].grid(True, linestyle="--")
for tick in axs[1].get_xticklabels() + axs[1].get_yticklabels():
    tick.set_fontweight('bold')
leg = axs[1].legend(title="Protein")
leg.get_title().set_fontweight('bold')

# --- SASA ---
sns.lineplot(data=avg_sasa, x="pH", y="sasa", hue="protein",
             marker="o", palette=custom_palette, ax=axs[2])
axs[2].set_title("Average SASA vs pH", fontweight='bold')
axs[2].set_ylabel("SASA (nm²)", fontweight='bold')
axs[2].set_xlabel("pH", fontweight='bold')
axs[2].grid(True, linestyle="--")
for tick in axs[2].get_xticklabels() + axs[2].get_yticklabels():
    tick.set_fontweight('bold')
leg = axs[2].legend(title="Protein")
leg.get_title().set_fontweight('bold')

# --- H-bonds ---
sns.lineplot(data=avg_hbonds, x="pH", y="avg_hbonds", hue="protein",
             marker="o", palette=custom_palette, ax=axs[3])
axs[3].set_title("Average H-bonds vs pH", fontweight='bold')
axs[3].set_ylabel("H-bonds", fontweight='bold')
axs[3].set_xlabel("pH", fontweight='bold')
axs[3].grid(True, linestyle="--")
for tick in axs[3].get_xticklabels() + axs[3].get_yticklabels():
    tick.set_fontweight('bold')
leg = axs[3].legend(title="Protein")
leg.get_title().set_fontweight('bold')

plt.tight_layout()

# ============================
# Save plots
# ============================
out_file_png = os.path.join(plot_folder, "combined_JAM_metrics_RMSD_first.png")
plt.savefig(out_file_png, dpi=300)

out_file_svg = os.path.join(plot_folder, "combined_JAM_metrics_RMSD_first.svg")
plt.savefig(out_file_svg)

plt.show()
print(f"✅ Plots saved to:\n- PNG: {out_file_png}\n- SVG: {out_file_svg}")

