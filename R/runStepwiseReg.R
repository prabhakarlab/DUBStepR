#' @title Run step-wise regression to order the features
#' @param ggc gene-gene correlation matrix
#' @return optimal feature set
#'
#'
runStepwiseReg <- function(ggc) {

    # Initialize variables
    ggc.cols <- ncol(ggc)
    ggc.rows <- nrow(ggc)
    ggc_centered <- ggc - Matrix::Matrix(data = Matrix::colMeans(ggc), ncol = ggc.cols, nrow = ggc.rows, byrow = TRUE)
    step = 1
    num_steps = 100
    num_regressions = 30
    step_seq = seq(from = step, to = num_steps, by = step)
    scree_values <- c(matrixcalc::frobenius.norm(x = as.matrix(ggc_centered)))
    names(scree_values) <- c("0")
    feature_genes <- c()
    
    # Progress bar
    message(" ")
    message("Running Stepwise Regression...")
    pb <- utils::txtProgressBar(style = 3)
    
    # For each step in the stepwise regression process
    for(i in step_seq) {
        if(i < num_regressions){
            # Set progress bar
            utils::setTxtProgressBar(pb = pb, value = i / num_regressions)
            
            # Compute GGC'*GGC
            #Crossprod saves time, but only at this step
            ggc_ggc <- Matrix::crossprod(x = ggc_centered, y = ggc_centered)
            dimnames(ggc_ggc) <- dimnames(ggc_centered)
            
            # Compute variance explained
            gcNormVec <- apply(ggc_ggc, 1, function(x) {
                matrixcalc::frobenius.norm(x)
            })
            
            # Compute norm of gene vectors
            gNormVec <- sqrt(abs(Matrix::diag(x = ggc_ggc)))
            
            # Select gene to regress out
            varExpVec <- gcNormVec / gNormVec
            regressed.genes <-
                names(sort(varExpVec, decreasing = TRUE))[1:step]
            
            # Add regressed gene to feature set
            feature_genes <- union(feature_genes, regressed.genes)
            
            # Obtain the variance explained by the regressed genes
            regressed.g <- ggc_centered[, regressed.genes]
            gtg <- as.numeric(t(regressed.g) %*% regressed.g)
            explained <-
                (regressed.g %*% t(ggc_ggc[regressed.genes, ])) / gtg
            eps <- ggc_centered - explained
            
            # Append variance explained value to vector for scree plot
            scree_values[paste0(i)] <- sum(varExpVec[regressed.genes])
            
            # Update the GGC to be the residual after regressing out these genes
            ggc_centered = eps

        } else {
            utils::setTxtProgressBar(pb = pb, value = i / num_regressions)
            scree_values[paste0(i)] <- scree_values[num_regressions]
            
        }
        
    }
    message(" ")
    message("Done.")
    
    # Plot scree plot to see change in Frobenius norm over genes
    scree_values <- scree_values[which(scree_values != 0)]
    
    # Find elbow point
    elbow_id <- findElbow(y = log(scree_values)[1:num_steps], ylab = "Log Variance Explained", plot = FALSE)
    
    # In case elbow is beyond set of feature genes
    if(elbow_id > length(feature_genes))
        elbow_id = length(feature_genes)
    
    # Initialise variables to add neighbours
    elbow_feature_genes <- feature_genes[1:elbow_id]
    neighbour_feature_genes <- elbow_feature_genes
    
    neighbour_fg_list <- as.list(elbow_feature_genes)
    names(neighbour_fg_list) <- elbow_feature_genes
    
    # Get list of GGC
    ggc_list <- apply(ggc, 2, function(x) {
        list(x)
    })
    ggc_list <- lapply(X = ggc_list, unlist)
    ggc_list <- lapply(X = ggc_list,
                       FUN = sort,
                       decreasing = TRUE)
    ggc_list <- lapply(ggc_list, function(x) {
        x[!(names(x) %in% neighbour_feature_genes)]
    })
    
    # Select potential candidates for next neighbour, changed to mclapply for faster computing
    candidateGGCList <-
        parallel::mclapply(ggc_list[neighbour_feature_genes], function(x) {
            x[which.max(x)]
        })
    
    # Progress bar
    message(" ")
    message("Adding correlated features...")
    max.ggc <- ggc.rows - length(neighbour_feature_genes)
    pb <- utils::txtProgressBar(min = 1, max = max.ggc, style = 3)
    
    # Adding neighbours of each gene
    for (i in 1:max.ggc) {
        
        # Set progress bar
        utils::setTxtProgressBar(pb = pb, value = i)
        
        # Select potential candidates for next neighbour
        candidateGGCNames <-
            unlist(lapply(candidateGGCList, names))
        
        candidateGGCVec <-
            unlist(candidateGGCList, use.names = FALSE)
        names(candidateGGCVec) <- candidateGGCNames
        
        # Select candidate with largest correlation to feature set
        nearest.index <- which.max(candidateGGCVec)
        whoseCandidate <- names(candidateGGCList[nearest.index])
        nearestNeighbour <- candidateGGCVec[nearest.index]
        
        # Add new neighbour to feature set, and remove from next neighbour candidate list
        neighbour_feature_genes <-
            append(neighbour_feature_genes, names(nearestNeighbour))
        
        # Update GGC list
        ggc_list <- lapply(ggc_list, function(x) {
            x[!(names(x) == names(nearestNeighbour))]
        })
        
        # Update list of nearest neighbour candidates
        candidateGGCList <-
            lapply(ggc_list[neighbour_feature_genes], function(x) {
                x[which.max(x)]
            })
    }
    message(" ")
    message("Done.")
    
    # Return
    return(list("feature.genes" = neighbour_feature_genes, "elbow.pt" = elbow_id))
}