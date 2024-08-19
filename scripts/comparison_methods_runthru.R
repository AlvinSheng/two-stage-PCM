# Running through the 4 competitor algorithms: 
# k-means global with PCA, k-means local with PCA, k-means global with whole curve, k-means local with whole curve



# PCM, psi set to 0

if (do_pcm0) {
  
  ptm <- proc.time()
  
  my_res_pcfs_psi0 <- met_gibbs_potts(y_array, M, W,
                                      n_burn_in = 10000, n_iter = 20000, thin = 1,
                                      keep_only = c(1, 31), potts_samps = potts_samps, true_cluster_labs = cluster_labs, finer_dim = finer_dim,
                                      label_unswitch = T, mu_inform = NULL, psi_fix = 0)
  
  potts_clustering_model_psi0_timing <- (proc.time() - ptm)[3]
  cat("minutes elapsed", potts_clustering_model_psi0_timing / 60)
  
  potts_clustering_model_psi0_ARI <- my_res_pcfs_psi0$mean.ARI
  
} else{
  potts_clustering_model_psi0_timing <- NA
  potts_clustering_model_psi0_ARI <- NA
}



# k-means global, PCA

y_array_da <- y_array_mean_impute(y_array)

fv_mat <- apply(y_array_da, 2, function(mat) c(t(mat)))

ptm <- proc.time()
kmeans_res <- kmeans(fv_mat, M)
kmeans_global_PCA_timing <- (proc.time() - ptm)[3]

if (!is.null(finer_dim)) { # if cluster_labs was generated at a finer resolution
  cluster_mat <- matrix(kmeans_res$cluster, nrow = num_subj, ncol = num_gp, byrow = T)
  kmeans_global_PCA_ARI <- adjustedRandIndex(cluster_labs_finer(cluster_mat, x_dim_coarse = S_dim$x_dim, y_dim_coarse = S_dim$y_dim, finer_dim = finer_dim), 
                                             cluster_labs)
} else {
  kmeans_global_PCA_ARI <- adjustedRandIndex(kmeans_res$cluster, t(cluster_labs)) # need to transpose cluster_labs, due to the ordering of subjects and grid regions
}



# k-means local, PCA

ptm <- proc.time()
kmeans_local_PCA_ARI <- kmeans_local_performance(fv_mat, num_subj, M, cluster_labs, S_dim, finer_dim)
kmeans_local_PCA_timing <- (proc.time() - ptm)[3]





# k-means global, whole curve

regime_smooth_mat <- do.call("rbind", lapply(regime_pcfs_smoothed$pcfs_smoothed, function(mat) t(mat)))
regime_smooth_mat <- scale(regime_smooth_mat)

if (calc_intens) { # incorporating intensities into the PCF data, as last index
  intens_mat <- do.call("rbind", local_intens_mats_scaled)
  regime_smooth_mat <- cbind(regime_smooth_mat, intens_mat)
}

pcf_da <- pcf_mat_mean_impute(regime_smooth_mat, num_subj, num_gp)

ptm <- proc.time()
kmeans_res <- kmeans(pcf_da, M)
kmeans_global_whole_curve_timing <- (proc.time() - ptm)[3]

if (!is.null(finer_dim)) { # if cluster_labs was generated at a finer resolution
  cluster_mat <- matrix(kmeans_res$cluster, nrow = num_subj, ncol = num_gp, byrow = T)
  kmeans_global_whole_curve_ARI <- adjustedRandIndex(cluster_labs_finer(cluster_mat, x_dim_coarse = S_dim$x_dim, y_dim_coarse = S_dim$y_dim, finer_dim = finer_dim), 
                                                     cluster_labs)
} else {
  kmeans_global_whole_curve_ARI <- adjustedRandIndex(kmeans_res$cluster, t(cluster_labs)) # need to transpose cluster_labs, due to the ordering of subjects and grid regions
}



# k-means local, whole curve

ptm <- proc.time()
kmeans_local_whole_curve_ARI <- kmeans_local_performance(pcf_da, num_subj, M, cluster_labs, S_dim, finer_dim)
kmeans_local_whole_curve_timing <- (proc.time() - ptm)[3]


