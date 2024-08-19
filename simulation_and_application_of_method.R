
library(here)
source(here("scripts/pkg_list.R"))



# General parameters 

# Fixed (True) Number of Clusters
M <- 3

# number of subjects
num_subj <- 50

# parameter for the Potts model
psi <- 1.29

# generate about lambda1 points per grid region
lambda1 <- 16
prob_mat <- cbind(c(0.5, 0.5), c(0.8, 0.2)) # ratio of type 1 and type 2 cells, given columnwise. 

# the string to pinpoint the folders to save the intermediary and final results in
factor_string <- paste0("M", M, "_N", num_subj, "_psi", psi, "_lambda1_", lambda1)

# for random number generation
seed <- M + num_subj + psi*100 + lambda1

# retrieving the potts_samps needed to interpolate the canonical stat values
x_dim_length <- 12
y_dim_length <- 10
res_string <- paste0("res", y_dim_length, "x", x_dim_length)

# You may see the warning about serializing the XGBoost model below.
# For the XGBoost package version used (version 1.7.5.1), this warning can be ignored.

potts_samps <- readRDS(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "potts_samps_xgb.rds")))



set.seed(609)



# Assume that all the spatial point patterns are rescaled onto a 
# specific grid. Even if spatial point patterns are of different
# sizes, data augmentation will take care of the missing grid regions 
x_dim <- seq(0, 12, length.out = 13)
y_dim <- seq(0, 10, length.out = 11)
S_dim <- list(x_dim = x_dim, y_dim = y_dim)

num_gp <- (length(S_dim$x_dim) - 1) * (length(S_dim$y_dim) - 1)

calc_intens <- TRUE

# Since all the spatial point patterns have been rescaled onto the same 
# window, the sequence of distances (r-sequence) can be the same. 

# rmax depends on the shortest unit-width S_dim uses
rmax <- shortside(Frame(owin(xrange = c(S_dim$x_dim[1], S_dim$x_dim[2]),
                             yrange = c(S_dim$y_dim[1], S_dim$y_dim[2]))))/2 
# assuming that the interesting pairwise interactions occur at distances less than half the grid region width

nr <- 513

r <- seq(0, rmax, length = nr)

cat("The r-sequence goes from 0 to", rmax, "(inclusive) with", nr, "r-values in total.")



# for the getNeighbors() function
source(here("scripts/getNeighbors.R"))

# Simulating the data, calculating local pcfs
source(here("scripts/sim_helper_fns.R"))

source(here("scripts", paste0("monte_carlo_sim_M", M, ".R")), local = T)
saveRDS(cluster_labs, file = here("intermediary_data", paste0("true_cluster_labs.rds")))
saveRDS(regime_spp, file = here("intermediary_data", paste0("regime_spp.rds")))

source(here("scripts/real_data_helper_fns.R"))

source(here("scripts/calculate_aggregated_LISAs.R"), local = T)
saveRDS(regime_pcfs_smoothed, file = here("intermediary_data", paste0("pcfs_smoothed.rds")))

# Applying the methods on the local pcfs 
source(here("scripts/clustering_algorithm_util.R"))
source(here("scripts/spatial_markov_model.R"))

# percentage of variance PCA should account for 
pca_cutoff <- 0.80

finer_dim <- NULL

label_unswitch <- T

source(here("scripts/methods_runthru.R"), local = T)
saveRDS(my_res_pcfs, file = here("MCMC_output", paste0("model_res.rds")))

do_pcm0 <- T
source(here("scripts/comparison_methods_runthru.R"), local = T)

ARI_vec <- c(potts_clustering_model_ARI, potts_clustering_model_psi0_ARI, kmeans_global_PCA_ARI, kmeans_local_PCA_ARI,
             kmeans_global_whole_curve_ARI, kmeans_local_whole_curve_ARI)
names(ARI_vec) <- c("potts_clustering_model", "potts_clustering_model_psi0", "kmeans_global_PCA", "kmeans_local_PCA", 
                    "kmeans_global_whole_curve", "kmeans_local_whole_curve")

timing_vec <- c(potts_clustering_model_timing, potts_clustering_model_psi0_timing, kmeans_global_PCA_timing, kmeans_local_PCA_timing,
                kmeans_global_whole_curve_timing, kmeans_local_whole_curve_timing)
names(timing_vec) <- c("potts_clustering_model", "potts_clustering_model_psi0", "kmeans_global_PCA", "kmeans_local_PCA", 
                       "kmeans_global_whole_curve", "kmeans_local_whole_curve")

list(ARIs = ARI_vec, 
     timings = timing_vec)


