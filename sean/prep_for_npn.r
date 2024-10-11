# Load the required package
library(data.table)

# Define the function
normalize_counts_to_dt <- function(counts) {
  # Convert counts to a data frame
  counts_df <- as.data.frame(counts)
  
  # Add row names as the first column
  counts_df <- cbind(Gene = rownames(counts_df), counts_df)
  
  # Scale values to the range 0 to 1
  scaled_counts_df <- as.data.frame(lapply(counts_df[-1], function(x) {
    (x - min(x)) / (max(x) - min(x))
  }))
  
  # Combine the scaled values with the gene names
  counts_df_scaled <- cbind(Gene = counts_df$Gene, scaled_counts_df)
  
  # Convert to a data table
  counts_dt <- as.data.table(counts_df_scaled)
  
  return(counts_dt)
}