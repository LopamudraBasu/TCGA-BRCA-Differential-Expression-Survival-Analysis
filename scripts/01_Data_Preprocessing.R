############################################################
# Project: TCGA-BRCA Differential Expression and Survival Analysis
# Script: 01_Data_Preprocessing.R
# Author: Lopamudra Basu
#
# Description:
# This script imports TCGA-BRCA RNA-seq expression,
# clinical and survival datasets and performs the
# initial data quality checks before downstream analyses.
############################################################

###############################
# Load Expression Data
###############################

expr <- read.delim(
  "TCGA.BRCA.sampleMap_HiSeqV2.gz",
  check.names = FALSE
)

###############################
# Load Clinical Data
###############################

pheno <- read.delim(
  "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix",
  check.names = FALSE
)

###############################
# Load Survival Data
###############################

survival <- read.delim(
  "survival_BRCA_survival.txt",
  check.names = FALSE
)

###############################
# Dataset Overview
###############################

cat("Expression matrix dimensions:\n")
dim(expr)

cat("Clinical data dimensions:\n")
dim(pheno)

cat("Survival data dimensions:\n")
dim(survival)

############################################################
# Prepare Clinical Sample Information
############################################################

# Select sample ID and sample type information
sample_info <- pheno[, c("sampleID", "sample_type")]

# Keep only Primary Tumor and Solid Tissue Normal samples
sample_info <- subset(
  sample_info,
  sample_type %in% c(
    "Primary Tumor",
    "Solid Tissue Normal"
  )
)

############################################################
# Match Expression and Clinical Samples
############################################################

# Identify common samples between expression and clinical datasets
common_samples <- intersect(
  colnames(expr)[-1],
  sample_info$sampleID
)

# Subset expression matrix using matched samples
expr_matrix <- expr[, c("sample", common_samples)]

# Reorder clinical information to match expression matrix columns
sample_info <- sample_info[
  match(common_samples, sample_info$sampleID),
]

############################################################
# Create Expression Matrix
############################################################

# Convert expression values to numeric matrix
expr_mat <- as.matrix(expr_matrix[, -1])

# Assign gene symbols as row names
rownames(expr_mat) <- expr_matrix$sample

############################################################
# Create Sample Groups
############################################################

group <- factor(
  sample_info$sample_type,
  levels = c(
    "Solid Tissue Normal",
    "Primary Tumor"
  )
)

############################################################
# Quality Control
############################################################

# Verify that sample order is identical
stopifnot(all(colnames(expr_mat) == sample_info$sampleID))

cat("Matched samples:", length(common_samples), "\n")
cat("Tumor samples:", sum(group == "Primary Tumor"), "\n")
cat("Normal samples:", sum(group == "Solid Tissue Normal"), "\n")