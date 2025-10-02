import csv
from pymol import cmd

# Step 3a: read RMSF values into a dictionary
rmsf_dict = {}
with open('/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMA/avg_files/pdb_figure/rmsf_values.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        resi = int(row['residue'])
        rmsf = float(row['rmsf'])
        rmsf_dict[resi] = rmsf

# Step 3b: initialize all B-factors to 0
cmd.alter('all', 'b=0.0')

# Step 3c: assign RMSF values to residues
for resi, rmsf in rmsf_dict.items():
    cmd.alter(f'resi {resi}', f'b={rmsf}')

# Step 3d: rebuild the structure to update visualization
cmd.rebuild()

