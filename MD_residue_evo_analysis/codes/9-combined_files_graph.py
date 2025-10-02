import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

# =========================================================
# 1️⃣ Paths
# =========================================================
input_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMC/combined_dynamics_sensitivity_consurf.csv"
out_dir   = "/content/drive/MyDrive/MD/JAM_pH/JAMC/res_graph_1"
os.makedirs(out_dir, exist_ok=True)  # create folder if not exists

# =========================================================
# 2️⃣ Load data
# =========================================================
df = pd.read_csv(input_csv)
print("Columns:", df.columns.tolist())
print(df.head())

# Expected columns:
# ['residue','SEQ','residue_pi','SCORE',
#  'avg_rmsf_all_ph','delta_rmsf','avg_sasa_all_ph','delta_sasa']

# =========================================================
# 3️⃣ Correlation analysis
# =========================================================
metrics = ['avg_rmsf_all_ph','delta_rmsf','avg_sasa_all_ph','delta_sasa']
corr_results = {}

for col in metrics:
    rho, p = spearmanr(df['residue_pi'], df[col])
    corr_results[col] = (rho, p)

print("\nSpearman correlations with residue pI:")
for m,(r,p) in corr_results.items():
    print(f"  {m:20s}: rho = {r:6.3f}, p = {p:6.3e}")

rho_rmsf, p_rmsf = spearmanr(df['SCORE'], df['delta_rmsf'])
rho_sasa, p_sasa = spearmanr(df['SCORE'], df['delta_sasa'])
print(f"\nConSurf SCORE vs ΔRMSF : rho={rho_rmsf:.3f}, p={p_rmsf:.3e}")
print(f"ConSurf SCORE vs ΔSASA : rho={rho_sasa:.3f}, p={p_sasa:.3e}")

# =========================================================
# 4️⃣ Scatter plots
# =========================================================
sns.set(style="whitegrid", context="talk")

# pI vs ΔRMSF
plt.figure(figsize=(8,6))
sns.scatterplot(
    data=df, x='residue_pi', y='delta_rmsf',
    hue='SCORE', size='avg_sasa_all_ph',
    palette='coolwarm', sizes=(20,200), edgecolor='black'
)
plt.title("Residue pI vs ΔRMSF\nColor=ConSurf SCORE, Size=avg SASA")
plt.xlabel("Residue pI")
plt.ylabel("ΔRMSF (across pH)")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "pI_vs_deltaRMSF.png"), dpi=300)
plt.show()

# Conservation vs ΔRMSF
plt.figure(figsize=(7,5))
sns.regplot(data=df, x='SCORE', y='delta_rmsf', scatter_kws={'s':40}, line_kws={'color':'red'})
plt.title("ConSurf SCORE vs ΔRMSF")
plt.xlabel("ConSurf SCORE (conservation)")
plt.ylabel("ΔRMSF")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "consurf_vs_deltaRMSF.png"), dpi=300)
plt.show()

# Conservation vs ΔSASA
plt.figure(figsize=(7,5))
sns.regplot(data=df, x='SCORE', y='delta_sasa', scatter_kws={'s':40}, line_kws={'color':'blue'})
plt.title("ConSurf SCORE vs ΔSASA")
plt.xlabel("ConSurf SCORE (conservation)")
plt.ylabel("ΔSASA")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "consurf_vs_deltaSASA.png"), dpi=300)
plt.show()

# =========================================================
# 5️⃣ Group analysis by residue type (acidic / basic / neutral)
# =========================================================
def classify_residue_pi(pi):
    if pi < 6:
        return "Acidic"
    elif pi > 8:
        return "Basic"
    else:
        return "Neutral"

df['pI_group'] = df['residue_pi'].apply(classify_residue_pi)

# Define the desired group order
group_order = ["Acidic", "Basic", "Neutral"]

group_stats = df.groupby('pI_group')[['delta_rmsf','delta_sasa']].mean().reindex(group_order).round(4)
print("\nMean dynamic changes by residue pI group:")
print(group_stats)

# Boxplot ΔRMSF by group
plt.figure(figsize=(6,5))
sns.boxplot(data=df, x='pI_group', y='delta_rmsf', order=group_order, palette='Set2')
plt.title("ΔRMSF by Residue pI Group")
plt.xlabel("Residue pI Group")
plt.ylabel("ΔRMSF")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "deltaRMSF_by_pI_group.png"), dpi=300)
plt.show()

# Boxplot ΔSASA by group
plt.figure(figsize=(6,5))
sns.boxplot(data=df, x='pI_group', y='delta_sasa', order=group_order, palette='Set3')
plt.title("ΔSASA by Residue pI Group")
plt.xlabel("Residue pI Group")
plt.ylabel("ΔSASA")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "deltaSASA_by_pI_group.png"), dpi=300)
plt.show()

# =========================================================
# 6️⃣ Identify key residues (high dynamics + high conservation)
# =========================================================
high_dyn_conserved = df[
    (df['delta_rmsf'] > 0.05) & (df['SCORE'] > 1.0)
].sort_values('delta_rmsf', ascending=False)

key_outfile = os.path.join(out_dir, "key_dynamic_conserved_residues.csv")
high_dyn_conserved.to_csv(key_outfile, index=False)
print(f"\nKey dynamic & conserved residues saved to:\n{key_outfile}")

# =========================================================
# 7️⃣ Export correlation summary
# =========================================================
corr_summary = pd.DataFrame([
    {"metric": m, "spearman_rho": r, "p_value": p}
    for m,(r,p) in corr_results.items()
])
corr_summary.to_csv(os.path.join(out_dir, "pI_correlations_summary.csv"), index=False)
print(f"Correlation summary saved to:\n{os.path.join(out_dir,'pI_correlations_summary.csv')}")

