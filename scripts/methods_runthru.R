
pca_summ <- summary(regime_pca)$importance

num_pc <- unname(which(pca_summ[3, ] >= pca_cutoff)[1])

pc_scores <- agg_pca[, 1:num_pc, drop = F]

# scaling the PC scores to have sd = 1, to be comparable with intensities
pc_scores <- scale(pc_scores)



start_subj <- 1
subj_range <- start_subj:(start_subj + num_subj - 1)

y_array <- preprocess_pc_scores(pc_scores = pc_scores, start_subj = start_subj, num_subj = num_subj, num_gp = num_gp)

if (calc_intens) {
  intens_array <- abind(local_intens_mats_scaled, along = 3)
  intens_array <- aperm(intens_array, c(3, 2, 1))
  y_array <- abind(y_array, intens_array, along = 2)
}



# assuming rook adjacency for the Potts layer
neighbors <- getNeighbors(matrix(1, (length(S_dim$y_dim) - 1), (length(S_dim$x_dim) - 1)), c(2, 2, 0, 0))



num_gp <- nrow(neighbors)

no_nbr_idx <- num_gp + 1

W <- matrix(0, num_gp, num_gp)

for (i in 1:num_gp) {
  
  nbr_idxs <- neighbors[i, neighbors[i, ] != no_nbr_idx]
  
  W[i, nbr_idxs] <- 1
  
}



ptm <- proc.time()

my_res_pcfs <- met_gibbs_potts(y_array, M, W,
                               n_burn_in = 10000, n_iter = 20000, thin = 1,
                               keep_only = c(1, 31), potts_samps = potts_samps, true_cluster_labs = cluster_labs, finer_dim = finer_dim,
                               label_unswitch = T, mu_inform = NULL, psi_fix = NULL)

potts_clustering_model_timing <- (proc.time() - ptm)[3]
cat("minutes elapsed", potts_clustering_model_timing / 60, "\n")



potts_clustering_model_ARI <- my_res_pcfs$mean.ARI

print(potts_clustering_model_ARI)

cat("\n")
