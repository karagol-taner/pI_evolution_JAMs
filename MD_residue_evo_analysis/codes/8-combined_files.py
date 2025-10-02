import pandas as pd

# -------------------------------
# File paths
# -------------------------------
rmsf_input_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMA/rmsf_per_residue.csv"
sasa_input_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMA/residue_sasa.csv"
consurf_csv    = "/content/JAMA_consurf_scores.csv"

output_rmsf_csv         = "/content/drive/MyDrive/MD/JAM_pH/JAMA/rmsf_avg_per_ph.csv"
output_sasa_csv         = "/content/drive/MyDrive/MD/JAM_pH/JAMA/sasa_avg_per_ph.csv"
output_combined_csv     = "/content/drive/MyDrive/MD/JAM_pH/JAMA/combined_per_ph_metrics.csv"
output_rmsf_change_csv  = "/content/drive/MyDrive/MD/JAM_pH/JAMA/rmsf_delta_per_residue.csv"
output_sasa_change_csv  = "/content/drive/MyDrive/MD/JAM_pH/JAMA/sasa_delta_per_residue.csv"
output_rmsf_avg_csv     = "/content/drive/MyDrive/MD/JAM_pH/JAMA/rmsf_avg_all_ph.csv"
output_sasa_avg_csv     = "/content/drive/MyDrive/MD/JAM_pH/JAMA/sasa_avg_all_ph.csv"
output_combined_change_csv = "/content/drive/MyDrive/MD/JAM_pH/JAMA/combined_dynamics_sensitivity_consurf.csv"

# -------------------------------
# Load CSVs
# -------------------------------
rmsf_df = pd.read_csv(rmsf_input_csv)
sasa_df = pd.read_csv(sasa_input_csv)
consurf_df = pd.read_csv(consurf_csv)

# Clean names and strip sequences
consurf_df.columns = consurf_df.columns.str.strip()
consurf_df['SEQ'] = consurf_df['SEQ'].astype(str).str.strip()
if 'SEQ' in rmsf_df.columns:
    rmsf_df['SEQ'] = rmsf_df['SEQ'].astype(str).str.strip()
if 'SEQ' in sasa_df.columns:
    sasa_df['SEQ'] = sasa_df['SEQ'].astype(str).str.strip()

# -------------------------------
# Residue pI dictionary
# -------------------------------
aa_pi = {
    'A': 6.00, 'R': 10.76, 'N': 5.41, 'D': 2.77, 'C': 5.07,
    'E': 3.22, 'Q': 5.65, 'G': 5.97, 'H': 7.95, 'I': 6.02,
    'L': 5.98, 'K': 9.74, 'M': 5.74, 'F': 5.48, 'P': 6.30,
    'S': 5.68, 'T': 5.60, 'W': 5.89, 'Y': 5.66, 'V': 5.96
}

# -------------------------------
# Shift residue numbering
# -------------------------------
rmsf_df['residue'] = rmsf_df['residue'] + 29
sasa_df['residue'] = sasa_df['residue'] + 29

# -------------------------------
# Aggregate RMSF and SASA
# -------------------------------
avg_rmsf_df = rmsf_df.groupby(['pH','residue'], as_index=False)['rmsf'].mean()
avg_rmsf_df.rename(columns={'rmsf':'avg_rmsf'}, inplace=True)
avg_rmsf_df.to_csv(output_rmsf_csv, index=False)

avg_sasa_df = sasa_df.groupby(['pH','residue'], as_index=False)['sasa'].mean()
avg_sasa_df.rename(columns={'sasa':'avg_sasa'}, inplace=True)
avg_sasa_df.to_csv(output_sasa_csv, index=False)

# -------------------------------
# Clean ConSurf POS safely
# -------------------------------
consurf_df = consurf_df[consurf_df['POS'].apply(lambda x: str(x).strip().isdigit())].copy()
consurf_df['POS'] = consurf_df['POS'].astype(int)

# Exclude ConSurf residues beyond simulation residues
max_res = max(avg_rmsf_df['residue'].max(), avg_sasa_df['residue'].max())
consurf_df = consurf_df[consurf_df['POS'] <= max_res]

# -------------------------------
# Merge RMSF + SASA + ConSurf + pI
# -------------------------------
merged_df = pd.merge(avg_rmsf_df, avg_sasa_df, on=['pH','residue'], how='outer')
merged_df = pd.merge(merged_df, consurf_df, left_on='residue', right_on='POS', how='left')
merged_df['SEQ'] = merged_df['SEQ'].astype(str).str.strip()
merged_df['residue_pi'] = merged_df['SEQ'].map(aa_pi)

merged_df = merged_df[['pH','residue','SEQ','residue_pi','SCORE','avg_rmsf','avg_sasa']]
merged_df.to_csv(output_combined_csv, index=False)
print(f"[✓] Combined per pH saved -> {output_combined_csv}")

# -------------------------------
# ΔRMSF and ΔSASA
# -------------------------------
rmsf_delta_df = avg_rmsf_df.groupby('residue')['avg_rmsf'].agg(['min','max']).reset_index()
rmsf_delta_df['delta_rmsf'] = rmsf_delta_df['max'] - rmsf_delta_df['min']
rmsf_delta_df.to_csv(output_rmsf_change_csv, index=False)

sasa_delta_df = avg_sasa_df.groupby('residue')['avg_sasa'].agg(['min','max']).reset_index()
sasa_delta_df['delta_sasa'] = sasa_delta_df['max'] - sasa_delta_df['min']
sasa_delta_df.to_csv(output_sasa_change_csv, index=False)

# -------------------------------
# Average across all pH
# -------------------------------
rmsf_avg_all_ph = avg_rmsf_df.groupby('residue')['avg_rmsf'].mean().reset_index()
rmsf_avg_all_ph.rename(columns={'avg_rmsf':'avg_rmsf_all_ph'}, inplace=True)
rmsf_avg_all_ph.to_csv(output_rmsf_avg_csv, index=False)

sasa_avg_all_ph = avg_sasa_df.groupby('residue')['avg_sasa'].mean().reset_index()
sasa_avg_all_ph.rename(columns={'avg_sasa':'avg_sasa_all_ph'}, inplace=True)
sasa_avg_all_ph.to_csv(output_sasa_avg_csv, index=False)

# -------------------------------
# Combine Δ and average data + SCORE
# -------------------------------
delta_summary_df = (
    rmsf_delta_df[['residue','delta_rmsf']]
    .merge(sasa_delta_df[['residue','delta_sasa']], on='residue')
    .merge(rmsf_avg_all_ph, on='residue')
    .merge(sasa_avg_all_ph, on='residue')
    .merge(consurf_df[['POS','SEQ','SCORE']], left_on='residue', right_on='POS', how='left')
)

delta_summary_df['SEQ'] = delta_summary_df['SEQ'].astype(str).str.strip()
delta_summary_df['residue_pi'] = delta_summary_df['SEQ'].map(aa_pi)

# ✅ Include SCORE in final export
delta_summary_df = delta_summary_df[['residue','SEQ','residue_pi','SCORE',
                                     'avg_rmsf_all_ph','delta_rmsf',
                                     'avg_sasa_all_ph','delta_sasa']]

delta_summary_df.to_csv(output_combined_change_csv, index=False)
print(f"[✓] ΔRMSF/ΔSASA summary saved -> {output_combined_change_csv}")


