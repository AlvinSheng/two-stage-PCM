
##### Step 1: computing PCFs

len_pcfs <- length(regime_spp)



regime_local_pcfs <- vector(mode = "list", length = len_pcfs)
local_intens_mats <- vector(mode = "list", length = len_pcfs)

# bw_vec <- rep(NA, len_pcfs)

ptm <- proc.time()

for (nn in 1:len_pcfs) {
  
  S <- regime_spp[[nn]]
  
  # # how the bandwidth is set in pcf.ppp according to the formula for the half-width of the Epanechnikov kernel
  # intens <- npoints(S) / spatstat.geom::area(S)
  # stoyan <- 0.15
  # h <- stoyan / sqrt(intens)
  # bw <- h / sqrt(5)
  
  loc_pcf_list <- local_window_calc_pcfs(S, S_dim, bw = NULL, r, min_cells = 2)
  regime_local_pcfs[[nn]] <- do.call("cbind", loc_pcf_list)
  
  if (calc_intens) {
    local_intens_mats[[nn]] <- local_window_calc_intens(S, S_dim) # local_intens_mats will be appended to y_array in methods_runthru.R
  }
  
}

# scaling the intensity values in local_intens_mats
intens_rbind <- do.call("rbind", local_intens_mats)
intens_means <- apply(intens_rbind, 2, function(vec) mean(vec, na.rm = T))
intens_sd <- apply(intens_rbind, 2, function(vec) sd(vec, na.rm = T))

local_intens_mats_scaled <- lapply(local_intens_mats, function(mat) sweep(mat, 2, intens_means, "-"))
local_intens_mats_scaled <- lapply(local_intens_mats_scaled, function(mat) sweep(mat, 2, intens_sd, "/"))

proc.time() - ptm





##### Step 2: Square rooting, smoothing, and trimming the pcfs

regime_pcfs_smoothed <- smooth_local_pcfs(r_grid = r, trim_rvals = c(0, rmax), fv_obj_list = regime_local_pcfs)

# The above object will be used for the non-PCA methods in comparison_methods_runthru.R





##### Step 3: Calculate FPCA
# will be used for the Potts clustering model and 2 of the comparison methods

# regime_local_pcfs and regime_pcfs_smoothed are organized
# so that the curves are columnwise, in the order of the grid regions.
regime_smooth_mat <- do.call("rbind", lapply(regime_pcfs_smoothed$pcfs_smoothed, function(mat) t(mat)))
# note that curves in regime_smooth_mat are ordered based on the grid region number, and then the other subjects are combined. 

regime_pca <- prcomp(regime_smooth_mat[complete.cases(regime_smooth_mat),])



agg_pca <- matrix(NA, nrow = nrow(regime_smooth_mat), ncol = ncol(regime_smooth_mat))

agg_pca[complete.cases(regime_smooth_mat)] <- regime_pca$x


