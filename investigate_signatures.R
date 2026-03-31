# investigate_signatures.R
# by Feng, Feng
# This script analyzes the mutation signatures structure in hpv H&N cancer samples.
# It reads Assignment_Solution_Activities.txt and generates PCA, UMAP, and Clustering plots.

# --- 1. Load / Install necessary packages ---
required_packages <- c("ggplot2", "umap", "pheatmap", "dplyr", "tidyr",
   "tibble", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse=", "))
    install.packages(new_packages, repos = "http://cran.us.r-project.org", quiet = TRUE)
}

library(ggplot2)
library(umap)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(here)
library(ggrepel)

# --- 2. Load and Prepare the Data ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  input_file <- args[1]
  message("Using input file from command line: ", input_file)
} else {
  input_file <- here("SNS101/vcf_data_assignment/Assignment_Solution/Activities", "Assignment_Solution_Activities.txt")
  message("No arguments provided. Using default input file: ", input_file)
}

if (!file.exists(input_file)) {
  stop("Input file '", input_file, "'' not found in the current directory.")
}

message("Reading ", input_file, "...")
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

samples <- data$Samples
activities <- as.matrix(data[, -1])
rownames(activities) <- samples

# Remove signatures (columns) that have exactly zero counts across all samples
activities_filtered <- activities[, colSums(activities) > 0]

# Normalize counts to proportions within each sample (compositional data)
proportions <- sweep(activities_filtered, 1, rowSums(activities_filtered), "/")

message("Data loaded. Total samples: ", nrow(proportions), ", Active signatures: ", ncol(proportions))

# --- 3. K-Means Clustering ---
message("Running K-Means Clustering...")
set.seed(42)

# Determine number of clusters (k). For small datasets, 3 is a good starting point.
# You can adjust this value based on your specific dataset.
k <- min(3, nrow(proportions) - 1)
kmeans_res <- kmeans(proportions, centers = k, nstart = 25)

# Extract cluster assignments
cluster_assignments <- as.factor(kmeans_res$cluster)
message("K-Means grouped samples into ", k, " clusters.")


# --- 4. Principal Component Analysis (PCA) ---
message("Running PCA...")

# prcomp centers variables by default, but shouldn't scale compositional data by variance
pca_result <- prcomp(proportions, center = TRUE, scale. = FALSE)

# Extract PC1 and PC2 and their variance explained
pca_variance <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df$Cluster <- cluster_assignments[rownames(pca_df)]

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "PCA of Mutation Signatures",
       subtitle = paste("Colored by K-Means Clusters (k =", k, ")"),
       x = sprintf("PC1 (%.1f%% Variance)", pca_variance[1]),
       y = sprintf("PC2 (%.1f%% Variance)", pca_variance[2])) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("PCA_Signatures.pdf", plot = p_pca, width = 8, height = 6)
message("Saved PCA_Signatures.pdf")


# --- 5. UMAP ---
message("Running UMAP...")
set.seed(42)  # For reproducible UMAP results

n_samples <- nrow(proportions)
umap_config <- umap.defaults
if(n_samples < 15) {
   umap_config$n_neighbors <- max(2, n_samples - 1)
}

umap_config$n_neighbors <- 4

umap_result <- umap(proportions, config = umap_config)
umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(umap_df)
umap_df$Cluster <- cluster_assignments[rownames(umap_df)]

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, show.legend = FALSE, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "UMAP of Mutation Signatures",
       subtitle = paste("Colored by K-Means Clusters (k =", k, ")"),
       x = "UMAP 1",
       y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave("UMAP_Signatures.pdf", plot = p_umap, width = 8, height = 6)
message("Saved UMAP_Signatures.pdf")


# --- 6. Hierarchical Clustering (Heatmap & Dendrogram) ---
message("Generating Hierarchical Clustering Heatmap...")

# We use pheatmap to do hierarchical clustering (Euclidean distance & complete linkage by default)
# and display the signature proportions.
# We transpose the matrix so samples are columns and signatures are rows.
heatmap_data <- t(proportions)

# Add k-means clusters as column annotations on the heatmap
annotation_col <- data.frame(KMeans_Cluster = cluster_assignments)
rownames(annotation_col) <- rownames(proportions)

pdf("Hierarchical_Clustering_Heatmap.pdf", width = 10, height = 8)
pheatmap(heatmap_data, 
         scale = "none",        # Values are already proportions
         cluster_cols = TRUE,   # Cluster Samples hierarchically
         cluster_rows = TRUE,   # Cluster Signatures
         annotation_col = annotation_col,
         main = "Hierarchical Clustering of Signature Proportions",
         color = colorRampPalette(c("white", "cornflowerblue", "firebrick3"))(100),
         border_color = "grey",
         angle_col = 45)
dev.off()
message("Saved Hierarchical_Clustering_Heatmap.pdf")

message("--- All analyses complete ---")
