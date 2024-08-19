# Script for operations related to running the Potts clustering model or the other clustering algorithms



#' Data pre-processing function
#' Assumes that the PC scores are ordered based on the grid region number within each subject,
#' after which other subjects are combined row-wise. 
#' @param pc_scores matrix of pc_scores, whose number of columns is the desired number of scores to use in clustering algorithm
#' @param start_subj first subject to extract from pc_scores
#' @param num_subj number of subjects to extract from pc_scores
#' @param num_gp number of grid regions in each subject
preprocess_pc_scores <- function(pc_scores, start_subj, num_subj, num_gp) {
  
  subj_range <- start_subj:(start_subj + num_subj - 1)
  
  num_scores <- ncol(pc_scores)
  
  mat_list <- vector(mode = "list", length = num_scores)
  
  for (j in 1:num_scores) {
    
    mat_list[[j]] <- matrix(pc_scores[(1 + (start_subj - 1) * num_gp):(num_gp * num_subj + (start_subj - 1) * num_gp), j],
                            nrow = num_subj, ncol = num_gp, byrow = T)
    
  }
  
  # appending it into an array
  y_array <- abind(mat_list, along = 3)
  
  # reordering the dimensions
  y_array <- aperm(y_array, perm = c(1, 3, 2))
  
  return(y_array)
  
}



# mean imputing missing values for the k-means clustering algorithm, 
# where the mean is taken over the principal component scores for a given dimension, within the given subject.
# Potts clustering model will use data augmentation

y_array_mean_impute <- function(y_array) {
  
  y_array_da <- y_array
  
  for (sn in 1:dim(y_array)[1]) {
    for (sj in 1:dim(y_array)[2]) {
      vec_isna <- is.na(y_array_da[sn, sj, ])
      y_array_da[sn, sj, ][vec_isna] <- mean(y_array_da[sn, sj, ], na.rm = T)
    }
  }
  
  return(y_array_da)
  
}



# mean imputing missing values for the k-means clustering algorithm that uses the entire pcf 
# (rather than dimensionality reduction via PCA), 
# where the mean is taken over the function values for a given r-value, within the given subject.
# Potts clustering model will use data augmentation
#' @param pcf_mat matrix of pcfs, binded row-wise by grid region first, then subject
#' @param num_subj number of subjects in pcf_mat
#' @param num_gp number of grid regions in each subject
pcf_mat_mean_impute <- function(pcf_mat, num_subj, num_gp) {
  
  pcf_da <- pcf_mat
  
  for (sn in 1:num_subj) {
    for (sj in 1:ncol(pcf_mat)) {
      vec <- pcf_da[(1 + (sn - 1) * num_gp):(sn * num_gp), sj]
      vec_isna <- is.na(vec)
      pcf_da[(1 + (sn - 1) * num_gp):(sn * num_gp), sj][vec_isna] <- mean(vec, na.rm = T)
    }
  }
  
  return(pcf_da)
  
}



# Return random folds for cross-validation
cv_folds <- function(num_subj, num_gp, num_folds) {
  
  test_fold_idxs <- matrix(NA, nrow = num_subj, ncol = num_gp)
  
  for (i in 1:num_subj) {
    
    test_fold_idxs[i, ] <- createFolds(1:num_gp, k = num_folds, list = F)
    
  }
  
  return(test_fold_idxs)
  
}



# Return the posterior predictions for out-of-sample grid regions, interlacing predictions from multiple folds into one data array
cv_preds <- function(test_fold_idxs, cv_models_list) {
  
  # make a list of posterior predictions from all models
  preds_list <- lapply(cv_models_list, function(model) model$mean.y.mat)
  
  oos_preds <- array(NA, dim = dim(preds_list[[1]]))
  
  for (k in 1:length(cv_models_list)) {
    
    for (sj in 1:dim(preds_list[[1]])[2]) {
      oos_preds[, sj, ][test_fold_idxs == k] <- preds_list[[k]][, sj, ][test_fold_idxs == k]
    }
    
  }
  
  return(oos_preds)
  
}



# Return the predicted pcfs (along with predicted intensities) constructed from the posterior predictions for out-of-sample grid regions' PC scores
cv_preds_pcf <- function(test_fold_idxs, cv_models_list, pca_cutoff, num_types, regime_pca_cv_folds, scale_cv_folds, r) {
  
  num_subj <- nrow(test_fold_idxs)
  num_gp <- ncol(test_fold_idxs)
  num_folds <- length(cv_models_list)
  nr <- length(regime_pca_cv_folds[[1]]$sdev)
  
  # make a list of posterior predictions from all models
  preds_list <- lapply(cv_models_list, function(model) model$mean.y.mat)
  
  pcfs_smoothed <- vector(mode = "list", length = num_subj)
  local_intens_mats <- vector(mode = "list", length = num_subj)
  
  for (sn in 1:num_subj) {
    
    pcfs <- matrix(NA, nrow = nr, ncol = num_gp)
    intens_mat <- matrix(NA, nrow = num_gp, ncol = num_types)
    
    for (k in 1:num_folds) {
      
      oos_idx <- which(test_fold_idxs[sn, ] == k)
      
      for (so in oos_idx) {
        
        pred_scores <- preds_list[[k]][sn, , so]
        
        # unscaling the pred_scores
        pred_scores <- pred_scores * scale_cv_folds[[k]][2, ]
        pred_scores <- pred_scores + scale_cv_folds[[k]][1, ]
        
        # assumes pred_scores is organized to contain num_pc PCs in the beginning, and intensities following
        pca_summ <- summary(regime_pca_cv_folds[[k]])$importance
        num_pc <- unname(which(pca_summ[3, ] >= pca_cutoff)[1])
        
        centered_pcf <- regime_pca_cv_folds[[k]]$rotation[, 1:num_pc] %*% as.matrix(pred_scores[1:num_pc])
        
        # un-center the pcf with the means stored in regime_pca_cv_folds
        pcfs[, so] <- centered_pcf + regime_pca_cv_folds[[k]]$center
        
        intens_mat[so, ] <- pred_scores[-c(1:num_pc)]
        
      }
      
    }
    
    pcfs_smoothed[[sn]] <- pcfs
    local_intens_mats[[sn]] <- intens_mat
    
  }
  
  regime_pcfs_smoothed <- list(pcfs_smoothed = pcfs_smoothed, r = r)
  
  return(list(regime_pcfs_smoothed = regime_pcfs_smoothed, local_intens_mats = local_intens_mats))
  
}



# transform the test observations according to the given PCA and scaling
pca_scale <- function(regime_pcfs_smoothed, local_intens_mats, pca_cutoff, regime_pca, given_scale) {
  
  r <- regime_pcfs_smoothed$r
  pcfs_list <- regime_pcfs_smoothed$pcfs_smoothed
  
  num_subj <- length(pcfs_list)
  num_gp <- nrow(local_intens_mats[[1]])
  
  nr <- length(r)
  num_types <- ncol(local_intens_mats[[1]])
  
  # organized to contain num_pc PCs in the beginning, and intensities following
  pca_summ <- summary(regime_pca)$importance
  num_pc <- unname(which(pca_summ[3, ] >= pca_cutoff)[1])
  
  trans_vals <- array(NA, dim = c(num_subj, num_pc + num_types, num_gp))
  
  for (sn in 1:num_subj) {
    
    pcfs_subj <- pcfs_list[[sn]]
    intens_subj <- local_intens_mats[[sn]]
    
    # do PCA
    pcfs_subj <- sweep(pcfs_subj, 1, regime_pca$center, "-") # centering the curves prior to PCA
    pc_scores <- t(regime_pca$rotation[, 1:num_pc]) %*% pcfs_subj # num_pc x num_gp
    
    # combine the quantities
    pc_scores <- rbind(pc_scores, t(intens_subj))
    
    # scaling the quantities
    pc_scores <- sweep(pc_scores, 1, given_scale[1, ], "-")
    pc_scores <- sweep(pc_scores, 1, given_scale[2, ], "/")
    
    # storing a matrix for all of the subject's grid regions
    trans_vals[sn, , ] <- pc_scores
    
  }
  
  return(trans_vals)
  
}



# transform the test observations according to the respective fold's PCA and scaling
cv_test_pca_scale <- function(test_fold_idxs, regime_pcfs_smoothed, local_intens_mats, pca_cutoff, regime_pca_cv_folds, scale_cv_folds) {
  
  r <- regime_pcfs_smoothed$r
  pcfs_list <- regime_pcfs_smoothed$pcfs_smoothed
  
  num_subj <- nrow(test_fold_idxs)
  num_gp <- ncol(test_fold_idxs)
  num_folds <- length(unique(c(test_fold_idxs)))
  nr <- length(r)
  num_types <- ncol(local_intens_mats[[1]])
  
  trans_vals <- vector(mode = "list", length = num_folds)
  
  for (k in 1:num_folds) {
    
    # assumes pred_scores is organized to contain num_pc PCs in the beginning, and intensities following
    pca_summ <- summary(regime_pca_cv_folds[[k]])$importance
    num_pc <- unname(which(pca_summ[3, ] >= pca_cutoff)[1])
    
    # transforming all grid regions (not just hold-outs) according to regime_pca_cv_folds and scale_cv_folds
    trans_vals[[k]] <- array(NA, dim = c(num_subj, num_pc + num_types, num_gp))
    
    for (sn in 1:num_subj) {
      
      pcfs_subj <- pcfs_list[[sn]]
      intens_subj <- local_intens_mats[[sn]]
      
      # do PCA
      pcfs_subj <- sweep(pcfs_subj, 1, regime_pca_cv_folds[[k]]$center, "-") # centering the curves prior to PCA
      pc_scores <- t(regime_pca_cv_folds[[k]]$rotation[, 1:num_pc]) %*% pcfs_subj # num_pc x num_gp
      
      # combine the quantities
      pc_scores <- rbind(pc_scores, t(intens_subj))
      
      # scaling the quantities
      pc_scores <- sweep(pc_scores, 1, scale_cv_folds[[k]][1, ], "-")
      pc_scores <- sweep(pc_scores, 1, scale_cv_folds[[k]][2, ], "/")
      
      # storing a matrix for all of the subject's grid regions
      trans_vals[[k]][sn, , ] <- pc_scores
      
    }
    
  }
  
  return(trans_vals)
  
}



# For the varying resolution simulation study, translate the cluster labels for the coarse resolution
# to the corresponding cluster labels for the finer resolution
#' @param coarse_cluster_labs assumes coarse_cluster_labs are of dimension N x number of coarse grid regions
#' @param x_dim_coarse breakpoints of grid regions to which coarse_cluster_labs were assigned, for x-axis
#' @param y_dim_coarse breakpoints of grid regions to which coarse_cluster_labs were assigned, for y-axis
#' @param finer_dim list of breakpoints (x_dim, y_dim) for the finer resolution 
#' @return corresponding cluster labels for the finer resolution
cluster_labs_finer <- function(coarse_cluster_labs, x_dim_coarse, y_dim_coarse, finer_dim) {
  
  # finding the corresponding coarse grid region for each fine grid region (column-major order)
  raster_grid_coarse <- owin(xrange = c(x_dim_coarse[1], tail(x_dim_coarse, 1)), yrange = c(y_dim_coarse[1], tail(y_dim_coarse, 1)), 
                             mask=matrix(TRUE, length(y_dim_coarse) - 1, length(x_dim_coarse) - 1))
  
  x_dim <- finer_dim$x_dim
  y_dim <- finer_dim$y_dim
  raster_grid <- owin(xrange = c(x_dim[1], tail(x_dim, 1)), yrange = c(y_dim[1], tail(y_dim, 1)), 
                      mask=matrix(TRUE, length(y_dim) - 1, length(x_dim) - 1))
  
  pc_spp <- pixelcentres(raster_grid)
  
  fine_gp_coords <- coords(pc_spp)
  corresp_coarse_gp <- do.call("cbind", nearest.raster.point(fine_gp_coords$x, fine_gp_coords$y, raster_grid_coarse))
  
  
  
  # finding the corresponding cluster labels for the finer resolution
  cluster_labs_fine <- matrix(NA, nrow(coarse_cluster_labs), nrow(corresp_coarse_gp))
  
  for (sn in 1:nrow(coarse_cluster_labs)) {
    
    for (gp_idx in 1:nrow(corresp_coarse_gp)) {
      
      coarse_gp_mat_idx <- corresp_coarse_gp[gp_idx, ]
      cluster_labs_fine[sn, gp_idx] <- coarse_cluster_labs[sn, mat_idx2vec(coarse_gp_mat_idx[1], coarse_gp_mat_idx[2], length(y_dim_coarse) - 1)]
      
    }
    
  }
  
  return(cluster_labs_fine)
  
}



#' performs k-means on each subject separately; returns the average ARI 
#' @param fv_mat matrix of function values or PC scores to apply k-means to. Assumes that the functions/scores are binded row-wise by grid region first, then subject
#' @param num_subj number of subjects
#' @param M number of clusters to input to k-means algorithm
#' @param true_cluster_labs true cluster labels to determine the clustering performance
#' @param finer_dim in case the true_cluster_labs were generated at a finer resolution
#' @param S_dim resolution that kmeans assumes
#' @return the performance in terms of ARI, averaged across subjects
kmeans_local_performance <- function(fv_mat, num_subj, M, true_cluster_labs, S_dim, finer_dim = NULL) {
  
  num_gp <- (length(S_dim$x_dim) - 1) * (length(S_dim$y_dim) - 1)
  
  ARI_vec <- rep(NA, num_subj)
  
  for (sn in 1:num_subj) {
    kmeans_res <- kmeans(fv_mat[(1 + (sn - 1) * num_gp):(sn * num_gp), ], M)
    
    if (!is.null(finer_dim)) { # if true_cluster_labs was generated at a finer resolution
      cluster_row <- matrix(kmeans_res$cluster, nrow = 1)
      ARI_vec[sn] <- adjustedRandIndex(cluster_labs_finer(cluster_row, x_dim_coarse = S_dim$x_dim, y_dim_coarse = S_dim$y_dim, finer_dim = finer_dim), 
                                       true_cluster_labs[sn, ])
    } else {
      ARI_vec[sn] <- adjustedRandIndex(kmeans_res$cluster, true_cluster_labs[sn, ])
    }
  }
  
  return(mean(ARI_vec))
  
}



# calculate the maximum number of matching adjacent pairs given the grid y_length by x_length
get_max_adj_pairs <- function(x_length, y_length) {
  
  # assuming rook adjacency for the Potts layer
  neighbors <- getNeighbors(matrix(1, y_length, x_length), c(2, 2, 0, 0))
  
  num_gp <- nrow(neighbors)
  
  return(sum(neighbors != num_gp + 1) / 2)
  
}



# Function to return smoothed object that can predict EPI value for given lambda (using constrained smoothing package scam)
# Assumes that the names of epi_curve are of the form "psi0", "psi0.01", ..., that correspond to the EPI values.
epi_fun <- function(epi_curve) {
  
  # extracting psi sequence
  psi_seq <- as.numeric(str_remove(names(epi_curve), "psi"))
  
  # returning scam model
  return(scam(epi_curve ~ s(psi_seq, bs="mpi")))
  
}



# Function to return smoothed object that can predict psi given the sum of equal pairs
psi_from_samps <- function(psi_samps) {
  
  # extracting psi sequence
  psi_seq <- as.numeric(str_remove(names(psi_samps[[1]]), "psi"))
  
  # extracting samps estimates
  samps <- psi_samps[[1]]
  
  # returning scam model
  return(scam(psi_seq ~ s(samps, bs="mpi", k = 30)))
  
}



# helper function to get sum_{i \sim j} I(C_i = C_j) statistic given the cluster labels
# for a single subject (cluster) and adjacency matrix W
pairwise_same_cluster_count <- function(cluster, W.triplet) {
  
  nbr_pairs <- cbind(W.triplet[, 2], W.triplet[, 1]) # row, col index
  # omitting the repeats
  nbr_pairs <- nbr_pairs[nbr_pairs[, 1] < nbr_pairs[, 2], ]
  
  same_clust_count <- sum(cluster[nbr_pairs[, 1]] == cluster[nbr_pairs[, 2]])
  
  return(same_clust_count)
  
}


