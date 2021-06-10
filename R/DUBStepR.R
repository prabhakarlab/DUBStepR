#' @author ranjanb
#' @title DUBStepR - Obtain a list of feature genes to characterise cell types
#' @param input.data input gene expression matrix (genes x cells)
#' @param min.cells minimum number of cells to filter genes out
#' @param species species to use for gene filtering: "human" (default), "mouse" and "rat"
#' @param optimise.features Determine optimal feature set using density index
#' @param k number of nearest neighbours. Default is 10
#' @param num.pcs number of principal components to represent sc data. Default is 20
#' @param error Acceptable error margin for kNN computation. Default is 0, but is set to 1 for large datasets
#' @return Returns optimal feature set
#'
#' @importClassesFrom Matrix dgCMatrix
#' @examples
#' res <- DUBStepR(input.data = pbmc_norm_small_data)
#' 
#' @export
#'
DUBStepR <- function(input.data, min.cells = 0.05*ncol(input.data), species = "human", optimise.features = TRUE, k = 10, num.pcs = 20, error = 0) {
    
    # Message
    message(" ")
    message("Running DUBStepR...")
    
    # Filter genes
    filt.data <- getFilteredData(data = input.data, min.cells = min.cells, species = species)
    
    # Compute gene-gene correlation matrix
    ggc.out <- getGGC(data = filt.data)
    
    # Run stepwise regression
    swreg.out <- runStepwiseReg(ggc = ggc.out$ggc)
    corr.info <- data.frame(feature.genes = swreg.out$feature.genes, corr.range = ggc.out$corr.range[swreg.out$feature.genes])
    
    if(optimise.features) {
        
        if(ncol(filt.data) > 10000) {
            error = 1
        }
        
        # Obtain optimal feature set
        density.out <- getOptimalFeatureSet(filt.data = filt.data, ordered.genes = swreg.out$feature.genes, elbow.pt = swreg.out$elbow.pt, k = k, num.pcs = num.pcs, error = error)
        
        dubStepR.out <- list("corr.info" = corr.info, "elbow.pt" = swreg.out$elbow.pt, "optimal.feature.genes" = density.out$optimal.feature.genes, "density.index" = density.out$density.index)
    } else {
        # Return all correlated genes
        dubStepR.out <- list("corr.info" = corr.info, "elbow.pt" = swreg.out$elbow.pt)
    }
    
    # Return DUBStepR output list
    return(dubStepR.out)
}
