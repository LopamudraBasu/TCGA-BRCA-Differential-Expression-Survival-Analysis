############################################################
# Project: TCGA-BRCA Differential Expression and Survival Analysis
# Script: 05_Heatmap_Visualization.R
# Author: Lopamudra Basu
#
# Description:
# Generate a heatmap of the top 50 differentially
# expressed genes identified by limma.
############################################################

###############################
# Load Required Packages
###############################

library(pheatmap)

###############################
# Select Top 50 DEGs
###############################

top50_genes <- rownames(sig_deg)[1:50]

###############################
# Prepare Heatmap Matrix
###############################

heatmap_data <- expr_mat[
  rownames(expr_mat) %in% top50_genes,
]

###############################
# Sample Annotation
###############################

annotation_col <- data.frame(
  SampleType = sample_info$sample_type
)

rownames(annotation_col) <- sample_info$sampleID

###############################
# Annotation Colors
###############################

annotation_colors <- list(
  SampleType = c(
    "Primary Tumor" = "#D73027",
    "Solid Tissue Normal" = "#4575B4"
  )
)

###############################
# Generate Heatmap
###############################

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = "complete",
  fontsize_row = 8,
  fontsize = 10,
  main = "Top 50 Differentially Expressed Genes"
)

###############################
# Save Heatmap
###############################

png(
  filename = "TCGA_BRCA_Heatmap_Top50_DEGs.png",
  width = 2200,
  height = 1800,
  res = 300
)

pheatmap(
  heatmap_data,
  scale = "row",
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_method = "complete",
  fontsize_row = 8,
  fontsize = 10,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

###############################
# Summary
###############################

cat("\n====================================\n")
cat("Heatmap Generated Successfully\n")
cat("====================================\n")
cat("Genes visualized:", nrow(heatmap_data), "\n")
cat("Samples:", ncol(heatmap_data), "\n")