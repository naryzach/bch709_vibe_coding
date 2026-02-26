#!/usr/bin/env Rscript
#
# Analyze Gasch et al. 2000 yeast stress response data
# Compute coefficient of variation (CV) and create heatmap of top 200 genes
#

library(data.table)
library(ggplot2)
library(pheatmap)
library(viridisLite)
library(scales)

# Read the Gasch 2000 data
message("Loading Gasch et al. 2000 stress response data...")
data <- fread("data/gasch2000.txt", sep = "\t", header = TRUE)

# Display column names for inspection
message("Column names:")
print(colnames(data))

# Extract gene IDs from UID column
gene_ids <- data$UID

# Skip non-numeric columns: UID, NAME, DESCRIPTION, GWEIGHT
# Select only numeric columns (stress conditions)
numeric_cols <- colnames(data)[sapply(data, is.numeric)]
expr_data <- data[, ..numeric_cols]

message(sprintf("Processing %d genes across %d stress conditions", 
                nrow(expr_data), ncol(expr_data)))

# Compute statistics for each gene
mean_expr <- rowMeans(expr_data, na.rm = TRUE)
sd_expr <- apply(expr_data, 1, sd, na.rm = TRUE)

# Compute CV = sd / abs(mean)
# Handle cases where mean is very close to zero
cv <- ifelse(abs(mean_expr) < 1e-10, NA, sd_expr / abs(mean_expr))

# Create data frame with results
cv_results <- data.table(
  gene_id = gene_ids,
  mean_expr = mean_expr,
  sd_expr = sd_expr,
  cv = cv
)

# Remove genes with all NA or zero variance
cv_results <- cv_results[!is.na(cv) & sd_expr > 0]

message(sprintf("After filtering: %d genes with non-zero variance", nrow(cv_results)))

# Select top 200 genes by CV (descending)
cv_top200 <- cv_results[order(-cv)][1:200]

# Round to 4 decimal places and save
cv_top200_rounded <- cv_top200[, .(
  gene_id,
  mean_expr = round(mean_expr, 4),
  sd_expr = round(sd_expr, 4),
  cv = round(cv, 4)
)]

# Save to TSV file
fwrite(cv_top200_rounded, 
       "results/yeast_stress_cv_top200.tsv", 
       sep = "\t")

message("Results saved to results/yeast_stress_cv_top200.tsv")

# Print top 10 to console
message("\nTop 10 genes by CV:")
print(cv_top200_rounded[1:10])

# Create heatmap data
# Select rows for top 200 genes and all expression columns
top200_gene_ids <- cv_top200$gene_id

# Create a data.table with gene IDs as a column for proper subsetting
expr_data_with_ids <- cbind(gene_id = gene_ids, expr_data)

# Filter to top 200 genes and reorder by CV
heatmap_data <- expr_data_with_ids[match(top200_gene_ids, expr_data_with_ids$gene_id), ]

# Set gene IDs as row names and remove the ID column
rownames(heatmap_data) <- heatmap_data$gene_id
heatmap_data$gene_id <- NULL

# Convert to matrix for pheatmap
heatmap_matrix <- as.matrix(heatmap_data)

message("\nGenerating heatmap...")

# Create heatmap with specified parameters
png("results/yeast_stress_cv_top200_heatmap.png",
    width = 1800, height = 1200, res = 200)

pheatmap(heatmap_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = viridis(256),
         main = "Yeast stress response, CV top200 (Gasch et al. 2000)",
         fontsize = 8,
         angle_col = 90,
         margins = c(10, 6))

dev.off()

message("Heatmap saved to results/yeast_stress_cv_top200_heatmap.png")

# Print summary statistics
message("\n=== Summary Statistics ===")
message(sprintf("Total genes analyzed: %d", nrow(cv_results)))
message("Top 200 genes selected")
message(sprintf("CV range (top 200): %.4f - %.4f", 
                min(cv_top200_rounded$cv), max(cv_top200_rounded$cv)))
message(sprintf("Mean expression range: %.4f - %.4f",
                min(cv_top200_rounded$mean_expr), max(cv_top200_rounded$mean_expr)))

message("\nScript completed successfully!")
