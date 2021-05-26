#' @title Filter the dataset by removing lowly expressed genes and mitochondrial, spike-in and ribosomal genes
#' @param data gene expression matrix
#' @param min.cells gene expression matrix
#' @return filtered gene-expression matrix
#'
#'
getFilteredData <- function(data, min.cells = 0.05*ncol(data)) {
    
    # Print data dimensions
    message(paste("Dimensions of input data:", dim(data)[[1]], "x", dim(data)[[2]]))
    
    # Filter data by gene expression
    filt.data <- data[Matrix::rowSums(data > 0) > min.cells, ]
    
    message(" ")
    message("Expression Filtering Done.")
    
    # Remove mitochondrial, spike-in and ribosomal genes as they do not assist in cell type separation
    filt.data <- filt.data[!grepl(pattern = "^MT-|^ERCC-|^RPS|^RPL", x = rownames(filt.data), ignore.case = TRUE), ]
    
    # Remove pseudogenes
    filt.data <- filt.data[!(rownames(filt.data) %in% pseudo_genes),]
    
    message("Mitochondrial, Ribosomal and Pseudo Genes Filtering Done.")
    
    
    # Print filtered data dimensions
    message(paste("Dimensions of filtered data: ",  dim(filt.data)[[1]], " x ", dim(filt.data)[[2]]))
    
    return(filt.data)
}
