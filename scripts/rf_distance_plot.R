# rf_distance_mds_plot.R
# This script reads an RF distance matrix, performs classical MDS,
# and plots the trees in 2D space with ggplot2. Each axis label includes the variance explained.

# Load required package
library(ggplot2)

# --- Read the RF Distance Matrix ---
# Assumes "rf_distance_matrix.rfdist" is a whitespace-separated file containing only numbers.
rf <- read.table("rf_distance_matrix.rfdist", sep = " ")

# Manually assign the tree names in the same order as the matrix rows/columns.
tree_names <- c("fasttree.treefile",
                "C20fixed.treefile",
                "basic.treefile",
                "qpfam_c60_g_tree_ufboot.treefile",
                "qpfam_c60_g_tree_ufboot_100.treefile",
                "GTR_c60_g_mwopt_tree_ufboot_100.treefile",
                "GTR_c60_g_tree_ufboot.treefile")
rownames(rf) <- tree_names
colnames(rf) <- tree_names

# --- Perform Classical MDS ---
# Use cmdscale with 'eig = TRUE' to retrieve eigenvalues for variance calculation.
mds_res <- cmdscale(as.matrix(rf), k = 2, eig = TRUE)
mds_coords <- mds_res$points
eig <- mds_res$eig

# Compute the variance explained by each dimension, using only positive eigenvalues.
pos_eig <- eig[eig > 0]
var_exp1 <- (eig[1] / sum(pos_eig)) * 100
var_exp2 <- (eig[2] / sum(pos_eig)) * 100

# Create a data frame for plotting.
mds_df <- data.frame(Dim1 = mds_coords[, 1],
                     Dim2 = mds_coords[, 2],
                     Tree = rownames(rf))

# --- Plot with ggplot2 ---
p <- ggplot(mds_df, aes(x = Dim1, y = Dim2, label = Tree)) +
    geom_point(size = 3, color = "steelblue") +
    geom_text(vjust = -0.5, hjust = 0.5, size = 4, color = "firebrick") +
    labs(title = "MDS Plot of RF Distances",
         x = paste0("MDS Dimension 1 (", round(var_exp1, 1), "%)"),
         y = paste0("MDS Dimension 2 (", round(var_exp2, 1), "%)")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(p)
