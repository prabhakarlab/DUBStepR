#' @title Compute the correlation range values for all genes in the gene-gene correlation matrix.
#' @param correlation_matrix gene-gene correlation matrix
#' @return list of p-values, adjusted p-values and correlation ranges for each gene
#'
#'
getCorrelationRange <- function(correlation_matrix) {

    # Sort each column of the correlation matrix
    sorted_corr = apply(correlation_matrix, 2, sort, decreasing = TRUE)

    # Select the second largest correlation for each gene
    max_corr = sorted_corr[3, ]

    # Select the smallest correlation for each gene
    min_corr = sorted_corr[nrow(sorted_corr), ]

    # Compute the correlation range
    diff_corr = max_corr-(0.75*min_corr)
    # names(diff_corr) <- names(rank_2nd_corr)

    # Sort genes in descending order of correlation range
    diff_corr <- sort(diff_corr, decreasing = T)

    # Return as list
    return(list("range"=diff_corr))
}
