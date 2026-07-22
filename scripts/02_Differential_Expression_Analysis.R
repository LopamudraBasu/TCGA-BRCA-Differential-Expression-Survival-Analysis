############################################################
# Project: TCGA-BRCA Differential Expression Analysis
# Script: 02_Differential_Expression_Analysis.R
# Author: Lopamudra Basu
#
# Description:
# Identify differentially expressed genes (DEGs) between
# Primary Tumor and Solid Tissue Normal samples using
# the limma package.
############################################################

###############################
# Load Required Packages
###############################

library(limma)
library(ggplot2)
library(ggrepel)

###############################
# Create Design Matrix
###############################

group <- factor(sample_info$sample_type)

design <- model.matrix(~0 + group)
colnames(design) <- c("Tumor", "Normal")

###############################
# Fit Linear Model
###############################

fit <- lmFit(expr_mat, design)

contrast.matrix <- makeContrasts(
  Tumor - Normal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

###############################
# Differential Expression Analysis
###############################

deg_results <- topTable(
  fit2,
  number = Inf,
  adjust.method = "BH"
)


###############################
# Filter Significant DEGs
###############################

sig_deg <- subset(
  deg_results,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

###############################
# Save DEG Results
###############################

write.csv(
  sig_deg,
  "TCGA_BRCA_DEGs.csv",
  row.names = TRUE
)

###############################
# Summary Statistics
###############################

cat("Total genes analysed:", nrow(deg_results), "\n")
cat("Significant DEGs:", nrow(sig_deg), "\n")


###############################
# Top 20 Differentially Expressed Genes
###############################

top20_deg <- head(
  sig_deg[
    order(sig_deg$adj.P.Val),
  ],
  20
)

write.csv(
  top20_deg,
  "TCGA_BRCA_top20_DEGs.csv",
  row.names = TRUE
)


############################################################
# Volcano Plot
############################################################


deg_results$Gene <- rownames(deg_results)

deg_results$Significant <- ifelse(
  deg_results$adj.P.Val < 0.05 &
    abs(deg_results$logFC) > 1,
  "DEG",
  "Non-DEG"
)

# Top 10 most significant genes
top10 <- head(
  deg_results[order(deg_results$adj.P.Val), ],
  10
)

volcano_plot <- ggplot(
  deg_results,
  aes(
    x = logFC,
    y = -log10(adj.P.Val)
  )
) +
  
  geom_point(
    aes(color = Significant),
    size = 1.3,
    alpha = 0.75
  ) +
  
  scale_color_manual(
    values = c(
      "DEG" = "#D73027",
      "Non-DEG" = "grey75"
    )
  ) +
  geom_vline(
    xintercept = c(-1,1),
    linetype = "dashed",
    colour = "grey50",
    linewidth = 0.4
  )+

geom_hline(
  yintercept = -log10(0.05),
  linetype = "dashed",
  colour = "grey50",
  linewidth = 0.4
 )+
  
  geom_text_repel(
    data = top10,
    aes(label = Gene),
    color = "black",
    size = 3.5,
    fontface = "bold",
    box.padding = 0.6,
    point.padding = 0.3,
    force = 3,
    segment.color = "grey40",
    max.overlaps = Inf
  ) +
  
  labs(
    title = "Differential Gene Expression in TCGA-BRCA",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P-value)
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    legend.title = element_blank(),
    legend.position = "right"
  )

volcano_plot

ggsave(
  "TCGA_BRCA_Volcano.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 600
)

############################################################
# End of Script
############################################################
