# original version, for running one set of parameters at a time, M < 6
library(here)
source(here("scripts/pkg_list.R"))

# local job script for making the all the KNN models stored in the potts_samps_xgb list, the input to spatial_markov_model.R

M <- 3
x_dim <- 12
y_dim <- 10

res_string <- paste0("res", y_dim, "x", x_dim)



# read in original canonical statistics
ordered_combined <- as.data.frame(read_delim(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "ordered_grid_canon.txt")), delim = " ", show_col_types = FALSE))
grid_pts_ordered <- ordered_combined[, 1:(M+1)]
canon_vals_ordered <- ordered_combined[, (M+2):(2*(M + 1))]
colnames(grid_pts_ordered) <- c(paste0("alpha", 1:M), "psi")



# Setting up the xgboost interpolation

xgboost_mod_list <- vector(mode = "list", length = M+1)
param <- list(booster = "gbtree", eval_metric = "rmse")
for (i in 1:(M+1)) {
  dtrain <- xgb.DMatrix(data = as.matrix(grid_pts_ordered), label = canon_vals_ordered[, i])
  xgboost_mod_list[[i]] <- xgb.train(param, dtrain, nrounds=100, max_depth = 12)
}



# input for spatial_markov_model.R
potts_samps_xgb <- list(xgboost_mod_list[[M+1]], xgboost_mod_list[1:M])



saveRDS(potts_samps_xgb, file = here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "potts_samps_xgb.rds")))




