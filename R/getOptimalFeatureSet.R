#' @title Determine the optimal feature set using Density Index (DI)
#' @param filt.data log-transformed filtered gene-expression matrix
#' @param ordered.genes genes ordered after stepwise regression
#' @param elbow.pt Elbow point to start determining optimal feature set
#' @param k number of nearest neighbours for CI computation
#' @param num.pcs number of principal components to represent sc data. Default is 20.
#' @param error Acceptable error margin for kNN computation. Default is 0, but is set to 1 for large datasets.
#' @return optimal set of feature genes
#'
#'
getOptimalFeatureSet <- function(filt.data, ordered.genes, elbow.pt = 25, k = 10, num.pcs = 20, error = 0) {

    # Initialise variables
    mean_knn_vec <- c()
    minNumGenes = ""
    numStepsUnchangedMin = 0
    
    # Progress bar
    message("Determining optimal feature set...")
    pb <- utils::txtProgressBar(min = elbow.pt, max = length(ordered.genes), style = 3)
    
    # For each neighbour
    for(num_genes in seq(from = elbow.pt, to = length(ordered.genes), by = 25)) {
        # Initialise number of genes
        neighbour_feature_genes <- ordered.genes[1:num_genes]
        
        # Run PCA on the feature data
        log.feature.data <-
            filt.data[neighbour_feature_genes, ]
        
        suppressWarnings({
            temp.seurat <- Seurat::CreateSeuratObject(counts = log.feature.data)
        
            Seurat::VariableFeatures(temp.seurat) <- neighbour_feature_genes
        
            temp.seurat <-
            Seurat::ScaleData(object = temp.seurat,
                              features = neighbour_feature_genes,
                              verbose = FALSE)
            temp.seurat <-
            Seurat::RunPCA(object = temp.seurat,
                           features = neighbour_feature_genes,
                           verbose = FALSE)})
        
        
        pca.data <- as.matrix(temp.seurat@reductions$pca@cell.embeddings[, 1:min(num.pcs, ncol(temp.seurat@reductions$pca@cell.embeddings))])
        rownames(pca.data) <- colnames(log.feature.data)
        
        # Compute k-NN distance
        system.time(
            my.knn <- RANN::nn2(
                data = pca.data,
                k = (k + 1),
                treetype = "kd",
                searchtype = "standard",
                eps = error
            )
        )
        
        nn.dists <- my.knn$nn.dists
        rownames(nn.dists) <- rownames(pca.data)
        
        # Remove first column as it consists of zeroes
        nn.dists <- nn.dists[,-1]
        
        # Calculate length scale to normalise distances
        sdVec <- stats::na.omit(temp.seurat@reductions$pca@stdev[1:num.pcs])
        length_scale <- sqrt(sum(sdVec ^ 2))
        
        # Scale k-NN distances by length scale
        mean_nn_dist <- mean(x = nn.dists)
        scaled_mean_nn_dist <- mean_nn_dist / length_scale
        names(scaled_mean_nn_dist) <- num_genes
        
        mean_knn_vec <-
            append(mean_knn_vec, scaled_mean_nn_dist)
        
        # Check if the minima has been updated
        if (which.min(mean_knn_vec) != minNumGenes) {
            minNumGenes = which.min(mean_knn_vec)
            numStepsUnchangedMin = 0
        } else {
            numStepsUnchangedMin = numStepsUnchangedMin + 1
        }
        
        # Set progress bar
        utils::setTxtProgressBar(pb = pb, value = num_genes)
    }
    
    # Determine optimal feature set
    optimal_feature_genes <- ordered.genes[1:as.numeric(names(minNumGenes))]
    
    message(" ")
    message("Done.")
    
    return(list("optimal.feature.genes" = optimal_feature_genes, "density.index" = mean_knn_vec))
}
