import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Data from the table
labels = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
JAMA_JAMB = [0.8394629, 0.874546, 0.8166244, 0.8590077, 0.8340313, 0.8833609, 0.4530124, 0.6523397, 0.7518182,
0.8429708, 0.6568363, 0.8765369, 0.8444415, 0.3569415, 0.7559558, 0.7480981, 0.7263794, 0.865871,
0.933395, 0.8498277]
JAMB_JAMC = [0.7211387, 0.9852315, 0.9333274, 0.9375555, 0.9529493, 0.9405178, 0.7280754, 0.9444954, 0.8984889,
0.8406604, 0.68905, 0.9449151, 0.9182613, 0.8193384, 0.8880582, 0.8740332, 0.8933799, 0.9352707,
0.9914435, 0.9547976]
JAMA_JAMC = [0.6611426, 0.8614692, 0.7450934, 0.818111, 0.834056, 0.8689652, 0.4174055, 0.6502374, 0.7140847,
0.7959727, 0.9163494, 0.8455381, 0.8333685, 0.3748823, 0.7559558, 0.7364964, 0.6669663, 0.8226961,
0.9213462, 0.7996234]

# Outlier detection thresholds
q1_1_5iqr_values = [0.565585775, 0.747324513, 0.500372963]
q3_1_5iqr_values = [1.037806175, 1.062965813, 1.038858663]

# Define color scheme
colors = ["#a6cee3", "#1f78b4", "#b2df8a"]

# Function to determine if a value is an outlier
def is_outlier(value, q1_1_5iqr, q3_1_5iqr):
    return value < q1_1_5iqr or value > q3_1_5iqr

# Plotting bar graph with less spacing between amino acid types
plt.figure(figsize=(16, 7))
x = np.arange(len(labels)) * 1.85  # Reduced spacing between groups (was 2)
width = 0.4  # Increased bar width (was 0.4)

# Bar plots with outliers using markers
for i in range(len(labels)):
    # Plot bars
    plt.bar(x[i] - width, JAMA_JAMB[i], width=width, color=colors[0], edgecolor='black', linewidth=0.8)
    plt.bar(x[i], JAMB_JAMC[i], width=width, color=colors[1], edgecolor='black', linewidth=0.8)
    plt.bar(x[i] + width, JAMA_JAMC[i], width=width, color=colors[2], edgecolor='black', linewidth=0.8)

    # Add markers for outliers
    if is_outlier(JAMA_JAMB[i], q1_1_5iqr_values[0], q3_1_5iqr_values[0]):
        plt.scatter(x[i] - width, JAMA_JAMB[i] + 0.02, color='red', marker='*', s=40, zorder=3)

    if is_outlier(JAMB_JAMC[i], q1_1_5iqr_values[1], q3_1_5iqr_values[1]):
        plt.scatter(x[i], JAMB_JAMC[i] + 0.02, color='red', marker='*', s=40, zorder=3)

    if is_outlier(JAMA_JAMC[i], q1_1_5iqr_values[2], q3_1_5iqr_values[2]):
        plt.scatter(x[i] + width, JAMA_JAMC[i] + 0.02, color='red', marker='*', s=40, zorder=3)

# Setting x-axis limits to reduce blank space on sides
plt.xlim(x[0] - width - 0.9, x[-1] + width + 0.9)  # Tighter x-axis limits

# Formatting for better visibility
plt.xticks(ticks=x, labels=labels, rotation=0, fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel("Values", fontsize=14)
plt.xlabel("Aminoacids", fontsize=14)
plt.title("Correlations between proteins based on aminoacids", fontsize=16, fontweight='bold')

# Creating a custom legend with outlier marker
legend_patches = [
    mpatches.Patch(facecolor=colors[0], edgecolor='black', label="JAMA-JAMB"),
    mpatches.Patch(facecolor=colors[1], edgecolor='black', label="JAMB-JAMC"),
    mpatches.Patch(facecolor=colors[2], edgecolor='black', label="JAMA-JAMC"),
    plt.Line2D([0], [0], marker='*', color='red', markersize=10, linestyle='None', label="Outliers"),
    mpatches.Patch(facecolor='red', alpha=0.1, label="Negative Correlation")
]

# Position the legend outside the plot area
plt.legend(handles=legend_patches, loc="upper left", bbox_to_anchor=(1, 1), fontsize=12, frameon=True, edgecolor='black')

# Enhanced Grid Visibility
plt.grid(axis='y', linestyle="--", alpha=0.6, linewidth=0.8)

# Add padding to the right side for the legend, but reduce other margins
plt.tight_layout()
plt.subplots_adjust(right=0.85, left=0.05)  # Reduced left margin

# Show the plot
plt.show()

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Define labels for amino acids
labels = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

# Updated values
JAMA_JAMB = [0.4074866, 0.2967109, 0.5061814, 0.4598549, 0.2345158, 0.39305, 0.2393702, 0.1530246, 0.3586483, 0.5302457,
             0.08762657, 0.443899, 0.3618325, 0.09367089, 0.3964115, 0.3175799, 0.3898953, 0.4792511, 0.3927921, 0.4837407]
JAMB_JAMC = [0.6977799, 0.9414373, 0.8439195, 0.7975012, 0.8454167, 0.7454731, 0.6653091, 0.9036145, 0.7834235, 0.5210554,
             0.2887, 0.7929035, 0.7246127, 0.7915991, 0.7643952, 0.7197083, 0.7985038, 0.783853, 0.9423867, 0.8696834]
JAMA_JAMC = [0.1481468, -0.001938557, -0.08244703, 0.07156181, 0.2347886, 0.2395728, 0.1433071, 0.1369713, 0.1332973, 0.2997216,
             0.8486584, 0.1097133, 0.2732185, 0.1539114, 0.1148177, 0.2562964, 0.05839834, 0.07269942, -0.08672027, -0.07525352]

# Outlier thresholds
q1_1_5iqr_values = [0.03410735, 0.5420258, -0.183299619]
q3_1_5iqr_values = [0.69615635, 1.0256546, 0.487555211]

# Define color scheme
colors = ["#a6cee3", "#1f78b4", "#b2df8a"]

# Function to determine if a value is an outlier
def is_outlier(value, q1_1_5iqr, q3_1_5iqr):
    return value < q1_1_5iqr or value > q3_1_5iqr

# Create the plot with white background
plt.figure(figsize=(16, 7))
ax = plt.gca()
ax.set_facecolor('white')

x = np.arange(len(labels)) * 1.85  # Adjusted spacing
width = 0.4  # Bar width

# Add a horizontal line at y=0
plt.axhline(y=0, color='black', linestyle='-', linewidth=1, zorder=1)

# Plot bars and highlight outliers
for i in range(len(labels)):
    # Draw bars with bottom at 0 for proper grounding
    plt.bar(x[i] - width, JAMA_JAMB[i], width=width, color=colors[0], edgecolor='black', linewidth=0.8, zorder=2)
    plt.bar(x[i], JAMB_JAMC[i], width=width, color=colors[1], edgecolor='black', linewidth=0.8, zorder=2)
    plt.bar(x[i] + width, JAMA_JAMC[i], width=width, color=colors[2], edgecolor='black', linewidth=0.8, zorder=2)

    # Add outlier markers
    if is_outlier(JAMA_JAMB[i], q1_1_5iqr_values[0], q3_1_5iqr_values[0]):
        marker_y = JAMA_JAMB[i] + 0.02 if JAMA_JAMB[i] > 0 else JAMA_JAMB[i] - 0.02
        plt.scatter(x[i] - width, marker_y, color='red', marker='*', s=40, zorder=3)

    if is_outlier(JAMB_JAMC[i], q1_1_5iqr_values[1], q3_1_5iqr_values[1]):
        marker_y = JAMB_JAMC[i] + 0.02 if JAMB_JAMC[i] > 0 else JAMB_JAMC[i] - 0.02
        plt.scatter(x[i], marker_y, color='red', marker='*', s=40, zorder=3)

    if is_outlier(JAMA_JAMC[i], q1_1_5iqr_values[2], q3_1_5iqr_values[2]):
        marker_y = JAMA_JAMC[i] + 0.02 if JAMA_JAMC[i] > 0 else JAMA_JAMC[i] - 0.02
        plt.scatter(x[i] + width, marker_y, color='red', marker='*', s=40, zorder=3)

# Find min and max values for proper y-axis limits
min_val = min(min(JAMA_JAMB), min(JAMB_JAMC), min(JAMA_JAMC))
max_val = max(max(JAMA_JAMB), max(JAMB_JAMC), max(JAMA_JAMC))
buffer = (max_val - min_val) * 0.1  # 10% buffer

# Adjust plot limits
plt.xlim(x[0] - width - 0.9, x[-1] + width + 0.9)
plt.ylim(min_val - 0.05, max_val + buffer)

# Format axes
plt.xticks(ticks=x, labels=labels, rotation=0, fontsize=12, fontweight='bold')
plt.yticks(fontsize=12)
plt.ylabel("Correlation Value", fontsize=14, fontweight='bold')
plt.xlabel("Amino Acids", fontsize=14, fontweight='bold')
plt.title("Partial Correlations between Proteins based on Amino Acids", fontsize=16, fontweight='bold')

# Add shaded region for negative values
plt.axhspan(min_val - 0.05, 0, alpha=0.1, color='red', zorder=0)


# Create legend
legend_patches = [
    mpatches.Patch(facecolor=colors[0], edgecolor='black', label="JAMA-JAMB"),
    mpatches.Patch(facecolor=colors[1], edgecolor='black', label="JAMB-JAMC"),
    mpatches.Patch(facecolor=colors[2], edgecolor='black', label="JAMA-JAMC"),
    plt.Line2D([0], [0], marker='*', color='red', markersize=10, linestyle='None', label="Outliers"),
    mpatches.Patch(facecolor='red', alpha=0.1, label="Negative Correlation"),

]

# Adjust legend position
plt.legend(handles=legend_patches, loc="upper left", bbox_to_anchor=(1, 1), fontsize=12, frameon=True, edgecolor='black')

# Improve grid visibility
plt.grid(axis='y', linestyle="--", alpha=0.6, linewidth=0.8)

# Optimize layout
plt.tight_layout()
plt.subplots_adjust(right=0.85, left=0.05)

# Show the plot
plt.show()
