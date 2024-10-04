quantile_normalization <- function(sce){
    counts = counts(sce)
    #determine the ranks of the counts
    counts_rank <- apply(counts,2,rank,ties.method="min")
    #sort the original from lowest to highest
    counts_sorted <- data.frame(apply(counts, 2, sort))
    # calculate the means
    counts_mean <- apply(counts_sorted, 1, mean)
    # set up a function for applying means to rank
    index_to_mean <- function(my_index, my_mean){
      return(my_mean[my_index])
    }
    #create the quantile counts
    counts_quantile_final <- apply(counts_rank, 2, index_to_mean, my_mean=counts_mean)
    rownames(counts_quantile_final) <- rownames(counts)
    #assigning the assay to the original object
    assay(sce, 4) <- counts_quantile_final
    return(sce)
}