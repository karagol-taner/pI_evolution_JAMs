import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Summary stats dict: (Count, Mean, StdDev)
data_summary = {
    'JAM2': {
        'All': (5643, 0.4487, 0.3164),
        'Other>Basic': (514, 0.4908, 0.3357),
        'Basic>Other': (720, 0.4170, 0.2596),
        'Other>Acidic': (536, 0.5133, 0.3304),
        'Acidic>Other': (522, 0.5156, 0.3095),
        'Basic->Acidic': (80, 0.4235, 0.2746),
        'Acidic->Basic': (58, 0.4191, 0.2944),
    },
    'JAM1': {
        'All': (5662, 0.5100, 0.3012),
        'Other>Basic': (536, 0.5533, 0.3053),
        'Basic>Other': (540, 0.4754, 0.2403),
        'Other>Acidic': (540, 0.5821, 0.2955),
        'Acidic>Other': (504, 0.5456, 0.2966),
        'Basic->Acidic': (60, 0.4636, 0.2339),
        'Acidic->Basic': (56, 0.4410, 0.2703),
    },
    'JAM3': {
        'All': (5871, 0.5310, 0.3405),
        'Other>Basic': (538, 0.5812, 0.3585),
        'Basic>Other': (720, 0.4688, 0.2999),
        'Other>Acidic': (540, 0.6050, 0.3442),
        'Acidic>Other': (702, 0.5245, 0.3315),
        'Basic->Acidic': (80, 0.4814, 0.3038),
        'Acidic->Basic': (78, 0.4574, 0.3318),
    }
}

# Mutation categories to keep consistent order
categories = ['Acidic->Basic', 'Basic->Acidic', 'Other>Basic', 'Other>Acidic', 'Acidic>Other', 'Basic>Other', 'All']

# Simulate data points for each protein and category
np.random.seed(42)  # for reproducibility

rows = []
for protein, cat_stats in data_summary.items():
    for cat in categories:
        count, mean, std = cat_stats[cat]
        # simulate pathogenicity scores with clipping to [0,1]
        scores = np.clip(np.random.normal(loc=mean, scale=std, size=count), 0, 1)
        for score in scores:
            rows.append({'protein': protein, 'mutation_category': cat, 'pathogenicity_score': score})

# Create DataFrame
df = pd.DataFrame(rows)

# Two-way ANOVA: effect of protein and mutation category on pathogenicity score
model = ols('pathogenicity_score ~ C(protein) + C(mutation_category) + C(protein):C(mutation_category)', data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)

import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt
import seaborn as sns

# Assuming your data is in a DataFrame `df` with columns:
# 'protein' (categorical), 'mutation_category' (categorical), and 'pathogenicity_score' (numeric)

# Fit the model again (if needed)
model = ols('pathogenicity_score ~ C(protein) * C(mutation_category)', data=df).fit()

# Tukey HSD for protein groups
tukey_protein = pairwise_tukeyhsd(df['pathogenicity_score'], df['protein'])
print(tukey_protein.summary())

# Tukey HSD for mutation_category groups
tukey_mut_cat = pairwise_tukeyhsd(df['pathogenicity_score'], df['mutation_category'])
print(tukey_mut_cat.summary())

from scipy.stats import ttest_ind, mannwhitneyu

# Compare Acidic→Basic vs all others (across all proteins)
acidic_to_basic = df[df['mutation_category'] == 'Acidic->Basic']['pathogenicity_score']
others = df[df['mutation_category'] != 'Acidic->Basic']['pathogenicity_score']

# Print group stats
print(f"\nAcidic→Basic (N={len(acidic_to_basic)}): Mean = {acidic_to_basic.mean():.4f}, Median = {acidic_to_basic.median():.4f}")
print(f"Other mutations (N={len(others)}): Mean = {others.mean():.4f}, Median = {others.median():.4f}")

# T-test (assumes normality)
t_stat, t_p = ttest_ind(acidic_to_basic, others, equal_var=False)
print(f"\nT-test: t = {t_stat:.3f}, p = {t_p:.3e}")

# Mann–Whitney U test (non-parametric)
u_stat, u_p = mannwhitneyu(acidic_to_basic, others, alternative='two-sided')
print(f"Mann–Whitney U: U = {u_stat}, p = {u_p:.3e}")

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu

# Function to calculate Cohen's d
def cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    s1, s2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    s_pooled = ((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2)
    return (np.mean(group1) - np.mean(group2)) / np.sqrt(s_pooled)

# Function to analyze direction vs. others
def analyze_direction_vs_other(df, direction, protein=None):
    if protein:
        df = df[df['protein'] == protein]

    target = df[df['mutation_category'] == direction]['pathogenicity_score']
    others = df[df['mutation_category'] != direction]['pathogenicity_score']

    # Statistics
    mean1, median1 = np.mean(target), np.median(target)
    mean2, median2 = np.mean(others), np.median(others)
    t_stat, t_p = ttest_ind(target, others, equal_var=False)
    u_stat, u_p = mannwhitneyu(target, others, alternative='two-sided')
    d = cohens_d(target, others)

    print(f"\n=== {direction} vs. Others {'(' + protein + ')' if protein else '(All Proteins Combined)'} ===")
    print(f"{direction} (N={len(target)}): Mean = {mean1:.4f}, Median = {median1:.4f}")
    print(f"Other Mutations (N={len(others)}): Mean = {mean2:.4f}, Median = {median2:.4f}")
    print(f"T-test: t = {t_stat:.3f}, p = {t_p:.3e}")
    print(f"Mann–Whitney U: U = {u_stat:.1f}, p = {u_p:.3e}")
    print(f"Cohen's d = {d:.3f}")

# Directions to test
directions = ['Acidic->Basic', 'Basic->Acidic']

# Run for combined and each protein
for direction in directions:
    analyze_direction_vs_other(df, direction)
    for prot in df['protein'].unique():
        analyze_direction_vs_other(df, direction, protein=prot)

# Define the custom color palette
custom_palette = {'JAM1': '#1b9e77', 'JAM2': '#d95f02', 'JAM3': '#7570b3'}

# Interaction plot: mean pathogenicity by mutation category and protein
plt.figure(figsize=(12,7))
sns.pointplot(data=df, x='mutation_category', y='pathogenicity_score', hue='protein',
              dodge=True, markers=['o', 's', 'D'], capsize=.1, errwidth=1, palette=custom_palette, hue_order=['JAM1', 'JAM2', 'JAM3'])
plt.title('AlphaMissense Pathogenicity by Mutation Category and Protein', fontsize=14, fontweight='bold')
plt.ylabel('Mean Pathogenicity Score', fontsize=12, fontweight='bold')
plt.xlabel('Mutation Category', fontsize=12, fontweight='bold')
plt.xticks(rotation=45, fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.grid(axis='y', linestyle="--", alpha=0.6, linewidth=0.8)
plt.tight_layout()

# Save the plot as SVG
plt.savefig("interaction_plot.svg", format='svg')

plt.show()

