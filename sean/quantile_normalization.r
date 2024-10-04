quantile_normalization <- function(sce){
  
  ordered_results <- apply(counts, 2, function(x) {
    # Get the order of the indices from highest to lowest
    ordered_indices <- order(x, decreasing = TRUE)
    # Return a list of ordered values and the original indices
    list(ordered_values = x[ordered_indices], original_indices = ordered_indices)
  })
  
  ordered_matrix <- matrix(
    unlist(lapply(ordered_results, function(x) x$ordered_values)),
    nrow = nrow(counts),
    byrow = FALSE
  )
  
  # Compute the average for each row
  row_averages <- rowMeans(ordered_matrix)
  
  # Replace each row's values with the row's average
  average_matrix <- matrix(rep(row_averages, ncol(ordered_matrix)),
                           nrow = nrow(ordered_matrix), byrow = FALSE)
  
  # Reconstruct original matrix based on saved indices
  reconstructed_matrix <- matrix(NA, nrow = nrow(counts), ncol = ncol(counts))
  
  for (i in 1:ncol(counts)) {
    original_indices <- ordered_results[[i]]$original_indices
    reconstructed_matrix[original_indices, i] <- average_matrix[, i]
  }
  
  colnames(reconstructed_matrix) <- colnames(counts)
  rownames(reconstructed_matrix) <- rownames(counts)
  
  assay(sce, 'quantile_norm') <- reconstructed_matrix             
  
  return(sce)
}