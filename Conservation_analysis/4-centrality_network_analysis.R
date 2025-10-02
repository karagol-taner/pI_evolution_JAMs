# Load required library
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)
residue_data <- data.frame(Aa, Ab, Ac, Ca, Cb, Cc, Da, Db, Dc, Ea, Eb, Ec, Fa, Fb, Fc, Ga, Gb, Gc, Ha, Hb, Hc, Ia, Ib, Ic, Ka, Kb, Kc, La, Lb, Lc, Ma, Mb, Mc, Na, Nb, Nc, Pa, Pb, Pc, Qa, Qb, Qc, Ra, Rb, Rc, Sa, Sb, Sc, Ta, Tb, Tc, Va, Vb, Vc, Wa, Wb, Wc, Ya, Yb, Yc)


# Compute correlation matrix for conservation scores
cor_matrix <- cor(residue_data, method = "pearson", use = "pairwise.complete.obs")  # Pearson correlation

# Lower correlation threshold to increase network connectivity
cor_threshold <- 0.1  # Allow more weak but meaningful connections
cor_matrix[abs(cor_matrix) < cor_threshold] <- 0  # Remove weak connections

# Convert correlation matrix into an adjacency matrix
adjacency_matrix <- abs(cor_matrix)  # Keep only strong correlations

# Convert the adjacency matrix into a graph object
network_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)


# Compute Betweenness Centrality for all residues
betweenness_scores <- betweenness(network_graph, directed = FALSE, weights = E(network_graph)$weight, normalized = TRUE)

# Compute Degree Centrality
degree_scores <- degree(network_graph, mode = "all", normalized = TRUE)

# Compute Weighted Degree Centrality (Strength)
strength_scores <- strength(network_graph, mode = "all")

# Check for network components (disconnected clusters)
components_info <- components(network_graph)
print(paste("Number of Disconnected Components:", components_info$no))

# Compute a score based on Betweenness Centrality (normalized to 0-100)
if (max(betweenness_scores) - min(betweenness_scores) > 0) {
  betweenness_score <- (betweenness_scores - min(betweenness_scores)) / (max(betweenness_scores) - min(betweenness_scores)) * 100
} else {
  betweenness_score <- betweenness_scores  # Avoid NaN if all values are the same
}

# Create a data frame of results
centrality_df <- data.frame(Residue = colnames(residue_data),
                            Betweenness_Centrality = betweenness_scores,
                            Degree_Centrality = degree_scores,
                            Weighted_Degree_Centrality = strength_scores,
                            Score = betweenness_score)

# Sort by highest Weighted Degree Centrality (instead of betweenness)
centrality_df <- centrality_df[order(-centrality_df$Weighted_Degree_Centrality), ]

# Print the ranked results
print(centrality_df)

# Save results to a CSV file
write.csv(centrality_df, "/content/centrality_scores.csv", row.names = FALSE)

# Plot network to visualize connectivity
plot(network_graph, vertex.size = 8, vertex.label.cex = 0.8, main = "Residue Interaction Network")

# Plot bar chart of residue importance based on Weighted Degree Centrality
barplot(centrality_df$Weighted_Degree_Centrality, names.arg = centrality_df$Residue, las = 2, col = "steelblue", main = "Residue Importance by Weighted Degree Centrality")


# Load necessary libraries
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("tidygraph", quietly = TRUE)) install.packages("tidygraph")

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)

# Compute the absolute adjacency matrix from the correlation matrix
adjacency_matrix <- abs(cor_matrix)  # Use your existing cor_matrix

# Convert to graph
network_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
E(network_graph)$weight <- E(network_graph)$weight * 5
tidy_graph <- as_tbl_graph(network_graph)

# ==== 1. Save to TIFF ====
svg("/content/residue_correlation_network.svg", width = 14, height = 14)

ggraph(tidy_graph, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), color = "black", width = 1, show.legend = FALSE) +
  geom_node_point(color = "#1f78b4", size = 7) +
  geom_node_text(aes(label = name), repel = TRUE, size = 5.5, fontface = "bold") +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1.2),  # Border around plot panel
    axis.text = element_text(face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    panel.grid.major = element_line(color = "grey70", linetype = "dotted", size = 0.6),
    panel.grid.minor = element_line(color = "grey70", linetype = "dotted", size = 0.6),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5)
  ) +
  ggtitle("Residue Correlation Network")

dev.off()

# Load required libraries
if (!requireNamespace("corrplot", quietly = TRUE)) install.packages("corrplot")

library(corrplot)  # For visualizing correlation matrix

# Assuming residue conservation scores are already stored in variables like Aa, Ab, Ac, etc.
# Create a data frame from residue conservation data
residue_data <- data.frame(Aa, Ab, Ac, Ca, Cb, Cc, Da, Db, Dc, Ea, Eb, Ec,
                           Fa, Fb, Fc, Ga, Gb, Gc, Ha, Hb, Hc, Ia, Ib, Ic,
                           Ka, Kb, Kc, La, Lb, Lc, Ma, Mb, Mc, Na, Nb, Nc,
                           Pa, Pb, Pc, Qa, Qb, Qc, Ra, Rb, Rc, Sa, Sb, Sc,
                           Ta, Tb, Tc, Va, Vb, Vc, Wa, Wb, Wc, Ya, Yb, Yc)

# Compute Pearson correlation matrix
cor_matrix <- cor(residue_data, method = "pearson", use = "pairwise.complete.obs")

# Apply correlation threshold (removes weak connections)
cor_threshold <- 0.1
cor_matrix[abs(cor_matrix) < cor_threshold] <- 0

# Print the correlation matrix in the console
print(cor_matrix)

svg("/content/correlation_matrix.svg", width = 14, height = 14)

corrplot(cor_matrix,
         method = "color",
         col = colorRampPalette(c("blue", "white", "red"))(200),
         tl.cex = 1.3,       # Increased text label size
         tl.col = "black",
         tl.font = 2,      # Set font to bold (2 for bold)
         tl.srt = 90,
         title = "",
         mar = c(1, 1, 3, 1),
         frame.plot = TRUE,
         addgrid.col = "gray",
         addshade = "all",
         shade.col = "gray90",
         cl.text.font = 2,     # Set font to bold for color legend text
         cl.cex = 1.3          # Increase size of color legend text
)


title(main = "Pearson Correlation Matrix", font.main = 2, cex.main = 2)

dev.off()
