#' @title Filter the dataset by removing lowly expressed genes and mitochondrial, spike-in and ribosomal genes
#' @param data gene expression matrix
#' @param min.cells gene expression matrix
#' @param species species to use for gene filtering: "human" (default), "mouse" and "rat"
#' @return filtered gene-expression matrix
#'
#'
getFilteredData <- function(data, min.cells = 0.05*ncol(data), species = "human") {
    
    # Print data dimensions
    message(paste("Dimensions of input data:", dim(data)[[1]], "x", dim(data)[[2]]))
    
    # Filter data by gene expression
    filt.data <- data[Matrix::rowSums(data > 0) > min.cells, ]
    
    message(" ")
    message("Expression Filtering Done.")
    
    # Remove mitochondrial, spike-in and ribosomal genes as they do not assist in cell type separation
    
    # Spike-in genes
    filt.data <- filt.data[!grepl(pattern = "^ERCC-", x = rownames(filt.data), ignore.case = TRUE), ]
    
    # Mitochondrial genes
    mito_table <- subset(allgene_list[[species]], grepl(pattern = "^MT-", x = allgene_list[[species]]$Gene.name, ignore.case = TRUE))
    mito_filter <- rownames(filt.data) %in% mito_table$Gene.stable.ID | rownames(filt.data) %in% mito_table$Gene.name | rownames(filt.data) %in% mito_table$Gene.Synonym
    
    # Ribosomal genes
    ribo_table <- subset(allgene_list[[species]], grepl(pattern = "^RPS|^RPL", x = allgene_list[[species]]$Gene.name, ignore.case = TRUE))
    ribo_filter <- rownames(filt.data) %in% ribo_table$Gene.stable.ID | rownames(filt.data) %in% ribo_table$Gene.name | rownames(filt.data) %in% ribo_table$Gene.Synonym
      
    filt.data <- filt.data[!(mito_filter | ribo_filter), ]
    
    # Remove pseudogenes
    ps_table <- subset(allgene_list[[species]], allgene_list[[species]]$Pseudogene.Status == "Pseudogene")
    
    ps_filter <- rownames(filt.data) %in% ps_table$Gene.stable.ID | rownames(filt.data) %in% ps_table$Gene.name | rownames(filt.data) %in% ps_table$Gene.Synonym
    filt.data <- filt.data[!ps_filter,]
    
    message("Mitochondrial, Ribosomal and Pseudo Genes Filtering Done.")
    
    
    # Print filtered data dimensions
    message(paste("Dimensions of filtered data: ",  dim(filt.data)[[1]], " x ", dim(filt.data)[[2]]))
    
    return(filt.data)
}
