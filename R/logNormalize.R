#' @title Log-transform and normalize data by sequencing depth
#' @param raw.data raw gene expression matrix
#' @param scale.factor scaling factor for normalization
#' @return log-normalized gene expression matrix
#'
logNormalize <- function(raw.data, scale.factor = 10000) {

    # Compute sequencing depth vector
    seqDepthVec <- Matrix::colSums(raw.data)
    
    # Scale sequencing depth by scaling factor
    seqDepthVec <- seqDepthVec/scale.factor
    
    # Normalise data by cell
    norm.data <- sapply(seq_along(seqDepthVec), function(index) {
        raw.data[,index]/seqDepthVec[index]
    })
    colnames(norm.data) <- colnames(raw.data)
    
    # Log-transform normalized data
    log.norm.data <- log(1+norm.data)
    
    message("Log-normalization Done.")
    
    # Return
    return(log.norm.data)
}
