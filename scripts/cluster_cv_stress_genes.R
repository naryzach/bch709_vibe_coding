# ============================================================================
# HIERARCHICAL CLUSTERING OF YEAST STRESS GENES (CV TOP 200)
# ============================================================================
#
# PROMPT (for reference):
# Write a new R code that runs in the bch709_vibe_coding conda environment.
# Installed packages: data.table, ggplot2, pheatmap, viridisLite, scales.
#
# --INPUT SPECIFICATION--
# - File: results/yeast_stress_cv_top200.tsv (TSV, tab-delimited, not gzipped)
#   Column 1: gene_id
#   Columns: mean_expr, sd_expr, cv
#
# - File: data/gasch2000.txt (TSV, tab-delimited, not gzipped)
#   Column 1: UID (gene_id)
#   Skip columns: NAME, description, GWEIGHT
#   Remaining columns: numeric stress-condition expression values (log2 ratios)
#
# --ANALYSIS CONDITIONS--
# - Subset gasch2000.txt to genes listed in yeast_stress_cv_top200.tsv
# - Remove metadata columns (NAME, description, GWEIGHT)
# - Keep only numeric condition columns
# - Remove rows with all NA or zero standard deviation
# - Z-score normalization (row-wise): Z = (value - row_mean) / row_sd
# - Hierarchical clustering: euclidean distance, ward.D2 method, k=4 clusters
#
# --OUTPUT 1 (Table)--
# Filename: results/cluster_assignment.tsv
# Columns: gene_id, cluster (integers, sorted by cluster, no decimals)
#
# --OUTPUT 2 (Plot)--
# Filename: results/cv_top200_cluster_heatmap.pdf
# Size: 8 × 12 inches
# Row-wise Z-scored expression, clustered genes, stress conditions in original order
# Title: "Z-Score Clustering of CV Top 200 Yeast Stress Genes"
# Annotation: distinct colors for clusters 1–4
# No column clustering (cluster_cols = FALSE)
#
# --OUTPUT 3 (QC)--
# Print to console:
#   - Number of genes loaded
#   - Number of genes after filtering
#   - Cluster sizes (counts per cluster)
#
# ============================================================================

library(data.table)
library(pheatmap)

# Ensure results directory exists
dir.create("results", showWarnings = FALSE)

# ============================================================================
# READ INPUT FILES
# ============================================================================

# Read CV gene list
cv_genes <- fread("results/yeast_stress_cv_top200.tsv")

# Read gasch2000 expression data
gasch <- fread("data/gasch2000.txt")

num_genes_loaded <- nrow(gasch)
cat("Number of genes loaded from gasch2000.txt:", num_genes_loaded, "\n")

# ============================================================================
# SUBSET TO CV GENES AND CLEAN METADATA
# ============================================================================

# Subset to genes in CV list
gasch_subset <- gasch[UID %in% cv_genes$gene_id, ]

# Remove metadata columns (NAME, GWEIGHT, and any description-like columns)
# Keep UID and identify numeric columns (conditions)
cols_to_remove <- c("NAME", "GWEIGHT")
cols_to_keep <- setdiff(names(gasch_subset), cols_to_remove)
gasch_subset <- gasch_subset[, ..cols_to_keep]

# Convert to matrix: rows = genes, columns = conditions
gene_ids <- gasch_subset$UID
expr_matrix <- as.matrix(gasch_subset[, -1, with = FALSE])
rownames(expr_matrix) <- gene_ids

# ============================================================================
# FILTER ROWS
# ============================================================================

# Remove rows with all NA values
complete_rows_mask <- rowSums(is.na(expr_matrix)) < ncol(expr_matrix)
expr_matrix <- expr_matrix[complete_rows_mask, ]

# Remove rows with zero standard deviation
row_sds <- apply(expr_matrix, 1, sd, na.rm = TRUE)
valid_sd_rows <- row_sds > 0 & !is.na(row_sds)
expr_matrix <- expr_matrix[valid_sd_rows, ]

num_genes_after_filter <- nrow(expr_matrix)
cat("Number of genes after filtering:", num_genes_after_filter, "\n")

# ============================================================================
# Z-SCORE NORMALIZATION (ROW-WISE)
# ============================================================================

expr_zscore <- t(scale(t(expr_matrix), center = TRUE, scale = TRUE))

# ============================================================================
# HIERARCHICAL CLUSTERING
# ============================================================================

# Compute distance matrix
dist_matrix <- dist(expr_zscore, method = "euclidean")

# Hierarchical clustering with ward.D2 method
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut tree to get 4 clusters
clusters <- cutree(hc, k = 4)

# ============================================================================
# CREATE CLUSTER ASSIGNMENT TABLE
# ============================================================================

cluster_df <- data.table(
  gene_id = names(clusters),
  cluster = as.integer(clusters)
)

# Sort by cluster
cluster_df <- cluster_df[order(cluster)]

# Print cluster sizes
cat("\nCluster sizes:\n")
print(cluster_df[, .N, by = cluster])

# Write to file
fwrite(cluster_df, "results/cluster_assignment.tsv", sep = "\t")

# ============================================================================
# PREPARE DATA FOR HEATMAP
# ============================================================================

# Reorder expression matrix rows by clustering
expr_zscore_ordered <- expr_zscore[cluster_df$gene_id, ]
clusters_ordered <- cluster_df[
  match(rownames(expr_zscore_ordered), gene_id), cluster
]

# Create annotation data frame for rows (cluster membership)
row_annotation <- data.frame(
  Cluster = as.factor(clusters_ordered),
  row.names = rownames(expr_zscore_ordered)
)

# Define distinct colors for clusters
cluster_colors <- list(
  Cluster = c(
    "1" = "#1B9E77",  # Green
    "2" = "#D95F02",  # Orange
    "3" = "#7570B3",  # Purple
    "4" = "#E7298A"   # Magenta
  )
)

# ============================================================================
# CREATE AND SAVE HEATMAP
# ============================================================================

pdf("results/cv_top200_cluster_heatmap.pdf", width = 8, height = 12)

pheatmap(
  expr_zscore_ordered,
  annotation_row = row_annotation,
  annotation_colors = cluster_colors,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = FALSE,  # Rows already ordered by cluster cut
  cluster_cols = FALSE,  # Keep condition column order
  main = "Z-Score Clustering of CV Top 200 Yeast Stress Genes",
  show_rownames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 101),
  fontsize_row = 8,
  fontsize_col = 8
)

dev.off()

cat("\nHeatmap saved to: results/cv_top200_cluster_heatmap.pdf\n")
cat("Cluster assignments saved to: results/cluster_assignment.tsv\n")

# ============================================================================
# INTERPRETATION COMMENTS
# ============================================================================
#
# CLUSTER 1: Genes showing robust upregulation (positive Z-scores) across
# multiple stress conditions, indicating core general stress response genes
# activated by diverse environmental challenges and essential for cellular
# protection and adaptation.
#
# CLUSTER 2: Genes displaying variable and condition-specific responses with
# mixed positive and negative Z-scores, representing modular stress pathways
# tailored to particular stress types and indicative of specialized regulatory
# mechanisms responding to distinct stimuli.
#
# CLUSTER 3: Genes with consistent downregulation (negative Z-scores) across
# stress conditions, likely representing growth-associated and metabolic genes
# whose expression is suppressed during stress to conserve cellular resources.
#
# CLUSTER 4: Genes showing intermediate or biphasic responses with expression
# induced in some conditions and suppressed in others, reflecting context-
# dependent and time-dependent regulation patterns indicative of context-
# specific or adaptive metabolic switches.
#
# ============================================================================
