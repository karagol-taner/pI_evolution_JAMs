import csv
from pymol import cmd

# Step 3a: read pi values into a dictionary
pi_dict = {}
with open('/Users/tanerkaragol/Desktop/JAM/MD/analysis/JAMC/avg_files/pdb_figure/pi_values.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        resi = int(row['residue'])
        pi = float(row['pi'])
        pi_dict[resi] = pi

# Step 3b: initialize all B-factors to 0
cmd.alter('all', 'b=0.0')

# Step 3c: assign pi values to residues
for resi, pi in pi_dict.items():
    cmd.alter(f'resi {resi}', f'b={pi}')

# Step 3d: rebuild the structure to update visualization
cmd.rebuild()

