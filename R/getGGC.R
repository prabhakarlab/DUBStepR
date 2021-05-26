#' @title Compute the correlation range values for all genes in the gene-gene correlation matrix
#' @param data log-transformed gene-expression matrix
#' @return list of genes with their z-transformed correlation range values
#'
#'
getGGC <- function(data) {
    message(" ")
    message("Computing GGC...")
    
    # Bin genes by mean expression
    num.bins = min(20, (nrow(data)-1))
    gene.mean <- Matrix::rowMeans(data)
    
    gene.bins <- cut(x = gene.mean, breaks = num.bins)
    names(gene.bins) <- names(gene.mean)
    
    # Compute correlation matrix
    t_data <- Matrix::t(data)
    correlation_matrix <- qlcMatrix::corSparse(X = t_data, Y = NULL)
    dimnames(correlation_matrix) <- list(rownames(data), rownames(data))
    
    message("Done.")
    
    # Obtain correlation range
    rangeObj <- getCorrelationRange(correlation_matrix = correlation_matrix)
    
    logRange <- log(1+rangeObj$range)
    
    # Z-transform correlation range
    zRangeList <- tapply(X = logRange, INDEX = gene.bins, FUN = function(x){
        (x - mean(x))/stats::sd(x)
    })
    
    if(all(is.na(zRangeList))) {
        stop("Feature correlation range could not be obtained.")
    }
    
    # Set NAs to 0
    zRangeList <- lapply(zRangeList, function(x){
        x[is.na(x)] <- 0
        return(x)
    })
    zRange <- unlist(zRangeList)
    names(zRange) <- unlist(lapply(strsplit(x = names(zRange), split = "\\]\\."), "[[", 2))
    
    # Obtain geneset with zRange > 0.7
    topGenes <- names(which(zRange > 0.7))
    
    # Reduce the GGC to use only these genes
    ggc <- correlation_matrix[topGenes, topGenes]
    
    # Return as matrix
    return(list("corr.range" = zRange, "ggc" = methods::as(ggc, "dgCMatrix")))
}
