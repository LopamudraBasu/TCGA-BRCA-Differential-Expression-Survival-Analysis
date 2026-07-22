############################################################
# Project: TCGA-BRCA Differential Expression and Survival Analysis
# Script: 04_Functional_Enrichment.R
# Author: Lopamudra Basu
#
# Description:
# Perform Gene Ontology (GO) and KEGG pathway enrichment
# analysis on significantly upregulated genes identified
# from TCGA-BRCA differential expression analysis.
############################################################

###############################
# Load Required Packages
###############################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

###############################
# Extract Upregulated Genes
###############################

up_genes <- rownames(
  subset(
    sig_deg,
    logFC > 1
  )
)

cat("Upregulated genes:", length(up_genes), "\n")

###############################
# Convert Gene Symbols to Entrez IDs
###############################

up_entrez <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat("Mapped genes:", nrow(up_entrez), "\n")

###############################
# GO Biological Process Enrichment
###############################

go_bp <- enrichGO(
  gene = up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

###############################
# Save GO Results
###############################

write.csv(
  as.data.frame(go_bp),
  "GO_BP_Enrichment.csv",
  row.names = FALSE
)

###############################
# GO Dot Plot
###############################

go_dotplot <- dotplot(
  go_bp,
  showCategory = 12,
  font.size = 11
) +
  ggtitle("GO Biological Process Enrichment") +
  theme_bw(base_size = 14)

go_dotplot

ggsave(
  "GO_Dot_Enrichment.png",
  plot = go_dotplot,
  width = 10,
  height = 9,
  dpi = 600
)

###############################
# GO Bar Plot
###############################

go_barplot <- barplot(
  go_bp,
  showCategory = 15,
  font.size = 11
)  +
  ggtitle("GO Biological Process Enrichment") +
  theme_bw(base_size = 14)

go_barplot

ggsave(
  "GO_Barplot.png",
  plot = go_barplot,
  width = 10,
  height = 9,
  dpi = 600
)

###############################
# KEGG Pathway Enrichment
###############################

kegg <- enrichKEGG(
  gene = up_entrez$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

###############################
# Convert Entrez IDs to Gene Symbols
###############################

kegg <- setReadable(
  kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

###############################
# Save KEGG Results
###############################

write.csv(
  as.data.frame(kegg),
  "KEGG_Enrichment.csv",
  row.names = FALSE
)

###############################
# KEGG Dot Plot
###############################

kegg_plot <- dotplot(
  kegg,
  showCategory = 12
) +
  ggtitle("KEGG Pathway Enrichment") +
  theme_bw(base_size = 14)

kegg_plot

ggsave(
  "KEGG_Enrichment.png",
  plot = kegg_plot,
  width = 10,
  height = 8,
  dpi = 600
)

###############################
# Summary
###############################

cat("\n====================================\n")
cat("Functional Enrichment Completed\n")
cat("====================================\n")

cat("Upregulated genes:", length(up_genes), "\n")
cat("Mapped genes:", nrow(up_entrez), "\n")
cat("GO terms:", nrow(go_bp), "\n")
cat("KEGG pathways:", nrow(kegg), "\n")