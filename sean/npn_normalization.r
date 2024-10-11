library('huge')

check_all_same <- function(x, my_tolerance = 1e-9){
  # This function returns TRUE if all the elements of the vector are the same
  # within a numerical tolerance levels
  # Thank you: https://stackoverflow.com/a/4752834
  #
  # Args:
  #   x: a numeric vector
  #   my_tolerance: how close must two numbers be for them to be considered equal?
  #
  # Returns:
  #   TRUE or FALSE
  if (is.numeric(x) & is.numeric(my_tolerance)) {
    return(all(abs(max(x) - min(x)) < my_tolerance))  
  } else {
    stop("Vector and tolerance given to check_all_same() must be numeric.")
  }
}

rescale_01 <- function(data_vector){
  # rescale values in a vector to [0,1]
  # Inputs: vector of numeric values
  # Returns: rescaled vector
  
  # if all the values are the same, return 0 vector
  if (check_all_same(data_vector)) {
    return(rep(0, length(data_vector)))
  } else {
    min_value <- min(data_vector)
    max_value <- max(data_vector)
    rescaled_values <- (data_vector - min_value)/(max_value - min_value)
    return(rescaled_values)
  }
}

rescale_datatable <- function(data_table, by_column = FALSE){
  # rescale each row (or column) of a data table to [0,1]
  # applies rescale_01() to each row (or column)
  # Inputs: gene expression data table
  #   first column of input is genes
  #   remaining columns are expression values
  #   by_column = FALSE rescale each row, if TRUE rescale each column
  # Returns: scaled gene expression data table
  
  data_table <- ensure_numeric_gex(data_table)
  
  data_matrix = data.matrix(data_table[, -1, with = F])
  
  # Rescale each row or column [0,1]
  if (by_column) {
    rescaled_data_matrix = apply(data_matrix, 2, rescale_01)
  } else {
    rescaled_data_matrix = t(apply(data_matrix, 1, rescale_01))
  }
  
  # Include gene symbols in result
  result = data.table(data.frame(data_table[,1], rescaled_data_matrix))
  colnames(result) <- colnames(data_table)
  
  result <- ensure_numeric_gex(result)
  
  return(data.table(result))
  
}

ensure_numeric_gex <- function(input_data) {
  # Numeric version of data.table
  #
  # Ensure gene expression values are numeric in a given data.table
  #
  # input_data: a data.table with gene in the first column and gene
  # expression in the remaining columns
  #
  # returns a data.table with numeric values for the gene expression columns
  
  if ("data.table" %in% class(input_data)) {
    
    # save gene names as character vector
    gene_names <- as.character(input_data[[1]])
    
    # force gene expression values to be numeric
    gex_values <- t(apply(input_data[,-1], 1, as.numeric))
    
    # create data frame of gene names and gene expression values
    # set column names to be same as input data
    return_df <- data.frame(gene_names, gex_values)
    names(return_df) <- names(input_data)
    
    # return as data.table
    return(data.table(return_df))
    
  } else {
    
    stop("\nInput must be a data.table")
    
  }
}


NPNSingleDT <- function(dt, zero.to.one = TRUE){
  # This function takes gene expression data in the form of a data.table
  # and returns a nonparanormal normalized, zero to one transformed
  # (if zero.to.one = TRUE) data.table
  #
  # Args:
  #   dt: data.table where the first column contains gene identifiers,
  #             the columns are samples, rows are gene measurements
  #	  zero.to.one: logical - should the data be zero to one transformed?
  #
  # Returns:
  #   npn.dt: nonparanormal normalized,
  #           if zero.to.one = TRUE zero to one transformed, data.table
  
  dt <- ensure_numeric_gex(dt)
  
  val <- data.frame(dt[, 2:ncol(dt), with = F])
  val.mat <- data.matrix(val)
  npn.mat <- huge.npn(t(val.mat), npn.func = "shrinkage",
                      npn.thresh = NULL,
                      verbose = FALSE)
  npn.dt <- data.table(cbind(as.character(dt[[1]]), t(npn.mat)))
  colnames(npn.dt) <- chartr(".", "-", colnames(dt))
  
  npn.dt <- ensure_numeric_gex(npn.dt)
  
  # npn.dt <- NAToZero(npn.dt)
  #  message("\tZero to one transformation...\n")
  if (zero.to.one) {
    npn.dt <- rescale_datatable(npn.dt)
  }
  return(data.table(npn.dt))
}