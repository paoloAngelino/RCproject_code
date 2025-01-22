# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)  # Get arguments passed to the script
p_thresh <- as.numeric(args[1])  # p-value threshold for differential expression
split_status <- tolower(args[2])  # Convert the second argument to lowercase for consistency

# Suppress warnings for cleaner output
options(warn = -1)

# Load required libraries
library(edgeR)
library(limma)
library(SingleCellExperiment)

# Define parameters
assay_to_use <- "scalelogcounts"  # Specify the assay to use for differential expression
file_path <- "/tmp/work/RCproject/GEO_singlecellexperiment.rds"  # Path to the SingleCellExperiment object

# Function to perform differential expression analysis using limma
run_limma <- function(spe, assay_name, verbose = FALSE, p.value = 0.05, lfc = 0.0) {
  # Create the design matrix
  design <- model.matrix(~0 + Response + batch, data = colData(spe))
  
  if (verbose) print(colnames(design))
  
  # Clean up column names in the design matrix
  colnames(design) <- gsub("^Response", "", colnames(design))
  colnames(design) <- gsub(" ", "_", colnames(design))
  
  if (verbose) print(colnames(design))
  
  # Define contrasts for differential expression analysis
  contr.matrix <- makeContrasts(
    Yes_vs_NO = yes - no,
    levels = colnames(design)
  )
  
  # Access assay data
  v <- assay(spe, assay_name)
  
  # Fit the linear model
  fit <- lmFit(v, design = design)
  fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)
  efit <- eBayes(fit_contrast, robust = TRUE)
  
  # Summarize results
  results_efit <- decideTests(efit, adjust.method = "none", p.value = p.value, lfc = lfc)
  summary_efit <- summary(results_efit)
  
  if (verbose) print(summary_efit)
  
  # Extract significant genes based on the p-value threshold
  my_de_results <- topTable(efit, coef = 1, sort.by = "none", n = Inf)
  sig_genes <- rownames(my_de_results[my_de_results$P.Value < p.value, ])
  
  return(sig_genes)
}

# Load the SingleCellExperiment object
sce <- readRDS(file_path)

# Subset the dataset based on the split_status argument
if (split_status == "true") {
  if (!file.exists("train_samples.txt")) {
    stop("train_samples.txt not found")
  }
  train_samples <- readLines("train_samples.txt")
  sce <- sce[, colnames(sce) %in% train_samples]
} else if (split_status == "false") {
  message("Running differential expression on the entire dataset")
} else {
  stop("Invalid value for split_status. Use 'true' or 'false'.")
}

# Subset data for different platforms/batches
all_data_sets <- sce
de4_data_sets <- sce[, colData(sce)$Platform %in% c("GPL13497", "GPL14951", "GPL15207", "GPL6102")]
rna_seq_data_sets <- sce[, colData(sce)$batch %in% c("GSE190826", "GSE209746")]

# Run differential expression analysis on the subsets
DE_alldatasets <- run_limma(all_data_sets, assay_name = assay_to_use, verbose = TRUE, p.value = p_thresh, lfc = 0.0)
DE_4datasets <- run_limma(de4_data_sets, assay_name = assay_to_use, verbose = TRUE, p.value = p_thresh, lfc = 0.0)
DE_bulkDatasets <- run_limma(rna_seq_data_sets, assay_name = assay_to_use, verbose = TRUE, p.value = p_thresh, lfc = 0.0)

# Find the intersection of results from all datasets
de_intersect_plus_bulk_genes <- intersect(union(DE_bulkDatasets, DE_4datasets), DE_alldatasets)

# Save the final gene set to an RDS file
saveRDS(de_intersect_plus_bulk_genes, "ann_gene_set.rds")
