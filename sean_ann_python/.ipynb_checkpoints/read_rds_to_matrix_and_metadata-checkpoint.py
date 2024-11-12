# Function to read RDS file and extract counts and metadata
def read_rds_to_matrix_and_metadata(file_path):
    # Load the RDS file in R
    ro.r(f"sce <- readRDS('{file_path}')")

    # Extract count data (assumed to be stored in assays)
    counts = ro.r('assay(sce, "scalelogcounts")')
    # Extract row (gene) and column (cell) names
    gene_names = ro.r('rownames(sce)')
    cell_names = ro.r('colnames(sce)')
    
    # Convert to a NumPy array
    counts_np = ro.conversion.rpy2py(counts)

    # Convert the counts matrix to a pandas DataFrame
    counts_df = pd.DataFrame(counts_np, index=gene_names, columns=cell_names)

    # Extract metadata from colData and convert to a pandas DataFrame directly
    metadata = ro.r('as.data.frame(colData(sce))')  # Get the colData as an R data frame
    metadata_df = pd.DataFrame(metadata)  # Convert R data frame to pandas DataFrame directly

    return counts_df, metadata_df