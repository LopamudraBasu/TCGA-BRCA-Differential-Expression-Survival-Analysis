############################################################
# Project: TCGA-BRCA Differential Expression and Survival Analysis
# Script: 03_Survival_Analysis.R
# Author: Lopamudra Basu
#
# Description:
# Perform Kaplan-Meier survival analysis for the top
# differentially expressed genes and identify prognostic
# biomarkers associated with overall survival.
############################################################

###############################
# Load Required Packages
###############################

library(survival)
library(survminer)

###############################
# Select Top 100 DEGs
###############################

top100_genes <- rownames(sig_deg)[1:100]

###############################
# Survival Analysis Function
###############################

survival_screen <- function(gene){
  
  # Extract gene expression
  gene_row <- expr[
    expr$sample == gene,
  ]
  
  # Skip if gene not found
  if(nrow(gene_row) == 0){
    return(
      data.frame(
        Gene = gene,
        Pvalue = NA
      )
    )
  }
  
  expression_df <- data.frame(
    sampleID = names(gene_row)[-1],
    expression = as.numeric(gene_row[1,-1])
  )
  
  # Merge with survival information
  survival_expression <- merge(
    expression_df,
    survival,
    by.x = "sampleID",
    by.y = "sample"
  )
  
  survival_expression <- survival_expression[
    !is.na(survival_expression$OS.time),
  ]
  
  # Skip if no samples remain
  if(nrow(survival_expression) == 0){
    return(
      data.frame(
        Gene = gene,
        Pvalue = NA
      )
    )
  }
  
  # High vs Low expression
  median_expression <- median(
    survival_expression$expression
  )
  
  survival_expression$group <- ifelse(
    survival_expression$expression >= median_expression,
    "High",
    "Low"
  )
  
  # Need two groups
  if(length(unique(survival_expression$group)) < 2){
    return(
      data.frame(
        Gene = gene,
        Pvalue = NA
      )
    )
  }
  
  fit <- survfit(
    Surv(OS.time, OS) ~ group,
    data = survival_expression
  )
  
  pvalue <- surv_pvalue(
    fit,
    data = survival_expression
  )$pval
  
  data.frame(
    Gene = gene,
    Pvalue = pvalue
  )
  
}

###############################
# Screen Top 100 Genes
###############################

survival_results <- do.call(
  rbind,
  lapply(
    top100_genes,
    survival_screen
  )
)

###############################
# Remove Missing Results
###############################

survival_results <- survival_results[
  !is.na(survival_results$Pvalue),
]

###############################
# Order by Significance
###############################

survival_results <- survival_results[
  order(survival_results$Pvalue),
]

###############################
# Significant Prognostic Genes
###############################

significant_genes <- subset(
  survival_results,
  Pvalue < 0.05
)

###############################
# Save Results
###############################

write.csv(
  survival_results,
  "TCGA_BRCA_Prognostic_Genes.csv",
  row.names = FALSE
)

###############################
# Summary
###############################

cat(
  "Genes screened:",
  nrow(survival_results),
  "\n"
)

cat(
  "Significant prognostic genes:",
  nrow(significant_genes),
  "\n"
)

print(significant_genes)

###############################
# Kaplan-Meier Plot
###############################

best_gene <- significant_genes$Gene[1]

cat(
  "Best prognostic gene:",
  best_gene,
  "\n"
)

gene_row <- expr[
  expr$sample == best_gene,
]

gene_expression <- data.frame(
  sampleID = names(gene_row)[-1],
  expression = as.numeric(gene_row[1,-1])
)

merged_gene <- merge(
  gene_expression,
  survival,
  by.x = "sampleID",
  by.y = "sample"
)

median_expression <- median(
  merged_gene$expression
)

merged_gene$group <- ifelse(
  merged_gene$expression >= median_expression,
  "High",
  "Low"
)

fit <- survfit(
  Surv(OS.time, OS) ~ group,
  data = merged_gene
)

km_plot <- ggsurvplot(
  fit,
  data = merged_gene,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#D73027","#2C7BB6"),
  title = paste(
    best_gene,
    "Expression and Overall Survival in TCGA-BRCA"
  ),
  xlab = "Time (days)",
  ylab = "Overall Survival Probability",
  legend.title = "",
  legend.labs = c(
    "High",
    "Low"
  )
)

print(km_plot)

ggsave(
  filename = paste0(
    best_gene,
    "_KaplanMeier.png"
  ),
  plot = km_plot$plot,
  width = 8,
  height = 6,
  dpi = 600
)

############################################################
# End of Script
############################################################