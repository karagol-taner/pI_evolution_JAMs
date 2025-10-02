import pandas as pd

print("=== JAM2 ===")

# Load the mutation file (adjust the filename and separator if needed)
df = pd.read_csv("JAM2.tsv", sep="\t")

# Descriptive statistics
mean_score = df['pathogenicity score'].mean()
std_dev = df['pathogenicity score'].std()
var = df['pathogenicity score'].var()
median = df['pathogenicity score'].median()
min_score = df['pathogenicity score'].min()
max_score = df['pathogenicity score'].max()
n = len(df)

# Output results
print("=== All mutations ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")

# Filter: FROM D or E → TO K or R (Acidic->Basic)
de_to_kr = df[
    (df['a.a.1'].isin(['D', 'E'])) &
    (df['a.a.2'].isin(['K', 'R']))
]

# Compute statistics
mean_score = de_to_kr['pathogenicity score'].mean()
std_dev = de_to_kr['pathogenicity score'].std()
median = de_to_kr['pathogenicity score'].median()
var = de_to_kr['pathogenicity score'].var()
min_score = de_to_kr['pathogenicity score'].min()
max_score = de_to_kr['pathogenicity score'].max()
n = len(de_to_kr)

# Print results
print("")
print("=== Acidic → Basic (D or E → K or R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter: FROM K or R → TO D or E (Basic->Acidic)
kr_to_de = df[
    (df['a.a.1'].isin(['K', 'R'])) &
    (df['a.a.2'].isin(['D', 'E']))
]

# Compute statistics
mean_score = kr_to_de['pathogenicity score'].mean()
std_dev = kr_to_de['pathogenicity score'].std()
median = kr_to_de['pathogenicity score'].median()
var = kr_to_de['pathogenicity score'].var()
min_score = kr_to_de['pathogenicity score'].min()
max_score = kr_to_de['pathogenicity score'].max()
n = len(kr_to_de)

# Print results
print("")
print("=== Basic → Acidic (K or R → D or E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")


# Filter: Other>Basic (FROM not K/R TO K or R)
basic_gain = df[
    ((df['a.a.2'] == 'K') | (df['a.a.2'] == 'R')) &      # mutated to K or R
    (~df['a.a.1'].isin(['K', 'R']))                      # original is not K or R
]

# Descriptive statistics
mean_score = basic_gain['pathogenicity score'].mean()
std_dev = basic_gain['pathogenicity score'].std()
var = basic_gain['pathogenicity score'].var()
median = basic_gain['pathogenicity score'].median()
min_score = basic_gain['pathogenicity score'].min()
max_score = basic_gain['pathogenicity score'].max()
n = len(basic_gain)

# Output results
print("")
print("=== Other>Basic Mutations (non-K/R → K or R only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter: Other>Acidic mutations
acidic_gain = df[
    ((df['a.a.2'] == 'D') | (df['a.a.2'] == 'E')) &      # mutated to D or E
    (~df['a.a.1'].isin(['D', 'E']))                      # original is not D or E
]

# Descriptive statistics
mean_score = acidic_gain['pathogenicity score'].mean()
std_dev = acidic_gain['pathogenicity score'].std()
var = acidic_gain['pathogenicity score'].var()
median = acidic_gain['pathogenicity score'].median()
min_score = acidic_gain['pathogenicity score'].min()
max_score = acidic_gain['pathogenicity score'].max()
n = len(acidic_gain)

# Output results
print("")
print("=== Other>Acidic Mutations (non-D/E → D or E only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter for Acidic>Other mutations: FROM D or E, TO anything EXCEPT D or E
acidic_loss = df[
    (df['a.a.1'].isin(['D', 'E'])) &                 # Original is D or E
    (~df['a.a.2'].isin(['D', 'E']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = acidic_loss['pathogenicity score'].mean()
std_dev = acidic_loss['pathogenicity score'].std()
median = acidic_loss['pathogenicity score'].median()
var = acidic_loss['pathogenicity score'].var()
min_score = acidic_loss['pathogenicity score'].min()
max_score = acidic_loss['pathogenicity score'].max()
n = len(acidic_loss)

# Print summary
print("")
print("=== Acidic>Other Mutations (D or E → non-D/E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter for Basic>Other mutations: FROM K or R, TO anything EXCEPT K or R
basic_loss = df[
    (df['a.a.1'].isin(['K', 'R'])) &                 # Original is Lysine or Arginine
    (~df['a.a.2'].isin(['K', 'R']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = basic_loss['pathogenicity score'].mean()
std_dev = basic_loss['pathogenicity score'].std()
median = basic_loss['pathogenicity score'].median()
var = basic_loss['pathogenicity score'].var()
min_score = basic_loss['pathogenicity score'].min()
max_score = basic_loss['pathogenicity score'].max()
n = len(basic_loss)

# Print summary
print("")
print("=== Basic>Other Mutations (K or R → non-K/R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

print("")

print("=== JAM1 ===")

# Load the mutation file (adjust the filename and separator if needed)
df = pd.read_csv("JAM1.tsv", sep="\t")

# Descriptive statistics
mean_score = df['pathogenicity score'].mean()
std_dev = df['pathogenicity score'].std()
var = df['pathogenicity score'].var()
median = df['pathogenicity score'].median()
min_score = df['pathogenicity score'].min()
max_score = df['pathogenicity score'].max()
n = len(df)

# Output results
print("=== All mutations ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")

# Filter: FROM D or E → TO K or R (Acidic->Basic)
de_to_kr = df[
    (df['a.a.1'].isin(['D', 'E'])) &
    (df['a.a.2'].isin(['K', 'R']))
]

# Compute statistics
mean_score = de_to_kr['pathogenicity score'].mean()
std_dev = de_to_kr['pathogenicity score'].std()
median = de_to_kr['pathogenicity score'].median()
var = de_to_kr['pathogenicity score'].var()
min_score = de_to_kr['pathogenicity score'].min()
max_score = de_to_kr['pathogenicity score'].max()
n = len(de_to_kr)

# Print results
print("")
print("=== Acidic → Basic (D or E → K or R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter: FROM K or R → TO D or E (Basic->Acidic)
kr_to_de = df[
    (df['a.a.1'].isin(['K', 'R'])) &
    (df['a.a.2'].isin(['D', 'E']))
]

# Compute statistics
mean_score = kr_to_de['pathogenicity score'].mean()
std_dev = kr_to_de['pathogenicity score'].std()
median = kr_to_de['pathogenicity score'].median()
var = kr_to_de['pathogenicity score'].var()
min_score = kr_to_de['pathogenicity score'].min()
max_score = kr_to_de['pathogenicity score'].max()
n = len(kr_to_de)

# Print results
print("")
print("=== Basic → Acidic (K or R → D or E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")


# Filter: Other>Basic (FROM not K/R TO K or R)
basic_gain = df[
    ((df['a.a.2'] == 'K') | (df['a.a.2'] == 'R')) &      # mutated to K or R
    (~df['a.a.1'].isin(['K', 'R']))                      # original is not K or R
]

# Descriptive statistics
mean_score = basic_gain['pathogenicity score'].mean()
std_dev = basic_gain['pathogenicity score'].std()
var = basic_gain['pathogenicity score'].var()
median = basic_gain['pathogenicity score'].median()
min_score = basic_gain['pathogenicity score'].min()
max_score = basic_gain['pathogenicity score'].max()
n = len(basic_gain)

# Output results
print("")
print("=== Other>Basic Mutations (non-K/R → K or R only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter: Other>Acidic mutations
acidic_gain = df[
    ((df['a.a.2'] == 'D') | (df['a.a.2'] == 'E')) &      # mutated to D or E
    (~df['a.a.1'].isin(['D', 'E']))                      # original is not D or E
]

# Descriptive statistics
mean_score = acidic_gain['pathogenicity score'].mean()
std_dev = acidic_gain['pathogenicity score'].std()
var = acidic_gain['pathogenicity score'].var()
median = acidic_gain['pathogenicity score'].median()
min_score = acidic_gain['pathogenicity score'].min()
max_score = acidic_gain['pathogenicity score'].max()
n = len(acidic_gain)

# Output results
print("")
print("=== Other>Acidic Mutations (non-D/E → D or E only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter for Acidic>Other mutations: FROM D or E, TO anything EXCEPT D or E
acidic_loss = df[
    (df['a.a.1'].isin(['D', 'E'])) &                 # Original is D or E
    (~df['a.a.2'].isin(['D', 'E']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = acidic_loss['pathogenicity score'].mean()
std_dev = acidic_loss['pathogenicity score'].std()
median = acidic_loss['pathogenicity score'].median()
var = acidic_loss['pathogenicity score'].var()
min_score = acidic_loss['pathogenicity score'].min()
max_score = acidic_loss['pathogenicity score'].max()
n = len(acidic_loss)

# Print summary
print("")
print("=== Acidic>Other Mutations (D or E → non-D/E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter for Basic>Other mutations: FROM K or R, TO anything EXCEPT K or R
basic_loss = df[
    (df['a.a.1'].isin(['K', 'R'])) &                 # Original is Lysine or Arginine
    (~df['a.a.2'].isin(['K', 'R']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = basic_loss['pathogenicity score'].mean()
std_dev = basic_loss['pathogenicity score'].std()
median = basic_loss['pathogenicity score'].median()
var = basic_loss['pathogenicity score'].var()
min_score = basic_loss['pathogenicity score'].min()
max_score = basic_loss['pathogenicity score'].max()
n = len(basic_loss)

# Print summary
print("")
print("=== Basic>Other Mutations (K or R → non-K/R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

print("")

print("=== JAM3 ===")

# Load the mutation file (adjust the filename and separator if needed)
df = pd.read_csv("JAM3.tsv", sep="\t")

# Descriptive statistics
mean_score = df['pathogenicity score'].mean()
std_dev = df['pathogenicity score'].std()
var = df['pathogenicity score'].var()
median = df['pathogenicity score'].median()
min_score = df['pathogenicity score'].min()
max_score = df['pathogenicity score'].max()
n = len(df)

# Output results
print("=== All mutations ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")

# Filter: FROM D or E → TO K or R (Acidic->Basic)
de_to_kr = df[
    (df['a.a.1'].isin(['D', 'E'])) &
    (df['a.a.2'].isin(['K', 'R']))
]

# Compute statistics
mean_score = de_to_kr['pathogenicity score'].mean()
std_dev = de_to_kr['pathogenicity score'].std()
median = de_to_kr['pathogenicity score'].median()
var = de_to_kr['pathogenicity score'].var()
min_score = de_to_kr['pathogenicity score'].min()
max_score = de_to_kr['pathogenicity score'].max()
n = len(de_to_kr)

# Print results
print("")
print("=== Acidic → Basic (D or E → K or R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter: FROM K or R → TO D or E (Basic->Acidic)
kr_to_de = df[
    (df['a.a.1'].isin(['K', 'R'])) &
    (df['a.a.2'].isin(['D', 'E']))
]

# Compute statistics
mean_score = kr_to_de['pathogenicity score'].mean()
std_dev = kr_to_de['pathogenicity score'].std()
median = kr_to_de['pathogenicity score'].median()
var = kr_to_de['pathogenicity score'].var()
min_score = kr_to_de['pathogenicity score'].min()
max_score = kr_to_de['pathogenicity score'].max()
n = len(kr_to_de)

# Print results
print("")
print("=== Basic → Acidic (K or R → D or E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")


# Filter: Other>Basic (FROM not K/R TO K or R)
basic_gain = df[
    ((df['a.a.2'] == 'K') | (df['a.a.2'] == 'R')) &      # mutated to K or R
    (~df['a.a.1'].isin(['K', 'R']))                      # original is not K or R
]

# Descriptive statistics
mean_score = basic_gain['pathogenicity score'].mean()
std_dev = basic_gain['pathogenicity score'].std()
var = basic_gain['pathogenicity score'].var()
median = basic_gain['pathogenicity score'].median()
min_score = basic_gain['pathogenicity score'].min()
max_score = basic_gain['pathogenicity score'].max()
n = len(basic_gain)

# Output results
print("")
print("=== Other>Basic Mutations (non-K/R → K or R only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter: Other>Acidic mutations
acidic_gain = df[
    ((df['a.a.2'] == 'D') | (df['a.a.2'] == 'E')) &      # mutated to D or E
    (~df['a.a.1'].isin(['D', 'E']))                      # original is not D or E
]

# Descriptive statistics
mean_score = acidic_gain['pathogenicity score'].mean()
std_dev = acidic_gain['pathogenicity score'].std()
var = acidic_gain['pathogenicity score'].var()
median = acidic_gain['pathogenicity score'].median()
min_score = acidic_gain['pathogenicity score'].min()
max_score = acidic_gain['pathogenicity score'].max()
n = len(acidic_gain)

# Output results
print("")
print("=== Other>Acidic Mutations (non-D/E → D or E only) ===")
print(f"Count: {n}")
print(f"Mean Pathogenicity Score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min: {min_score:.4f}, Max: {max_score:.4f}")

# Filter for Acidic>Other mutations: FROM D or E, TO anything EXCEPT D or E
acidic_loss = df[
    (df['a.a.1'].isin(['D', 'E'])) &                 # Original is D or E
    (~df['a.a.2'].isin(['D', 'E']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = acidic_loss['pathogenicity score'].mean()
std_dev = acidic_loss['pathogenicity score'].std()
median = acidic_loss['pathogenicity score'].median()
var = acidic_loss['pathogenicity score'].var()
min_score = acidic_loss['pathogenicity score'].min()
max_score = acidic_loss['pathogenicity score'].max()
n = len(acidic_loss)

# Print summary
print("")
print("=== Acidic>Other Mutations (D or E → non-D/E) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")

# Filter for Basic>Other mutations: FROM K or R, TO anything EXCEPT K or R
basic_loss = df[
    (df['a.a.1'].isin(['K', 'R'])) &                 # Original is Lysine or Arginine
    (~df['a.a.2'].isin(['K', 'R']))                  # Mutated TO something else
]

# Calculate statistics
mean_score = basic_loss['pathogenicity score'].mean()
std_dev = basic_loss['pathogenicity score'].std()
median = basic_loss['pathogenicity score'].median()
var = basic_loss['pathogenicity score'].var()
min_score = basic_loss['pathogenicity score'].min()
max_score = basic_loss['pathogenicity score'].max()
n = len(basic_loss)

# Print summary
print("")
print("=== Basic>Other Mutations (K or R → non-K/R) ===")
print(f"Count: {n}")
print(f"Mean pathogenicity score: {mean_score:.4f}")
print(f"Median: {median:.4f}")
print(f"Standard Deviation: {std_dev:.4f}")
print(f"Variance: {var:.4f}")
print(f"Min score: {min_score:.4f}")
print(f"Max score: {max_score:.4f}")


