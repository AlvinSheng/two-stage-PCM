# Script for creating the Monte Carlo data needed to compute the normalizing constant 
# in the spatial Potts model, to be used in spatial_markov_model.R, 
# now with offsets

# Focused on the psi parameter while accounting for alpha parameters

# The data necessary for computing the normalizing constant will be computed for 
# psi = 0.00, 0.01, 0.02,... , 2.50, for a given M for the corresponding alpha_2, ..., alpha_M, in the range [-5, 5] and the given grid x_dim x y_dim. 



# general code for the parameters passed in

set.seed(seed+1, "L", "I", "Re")

alpha_sd <- 4

# # version that uses the same alphas for every psi
# 
# grid_pre <- matrix(NA, nrow = nas, ncol = iclust)
# 
# for (i in 1:iclust) {
#   grid_pre[, i] <- rtruncnorm(n=nas, a=alpha_range[1], b=alpha_range[2], mean=0, sd=alpha_sd) # runif(n = nas, min = alpha_range[1], max = alpha_range[2])
# }
# 
# # including psi
# 
# grid_pts <- matrix(rep(grid_pre, each = length(psi_seq)), ncol = ncol(grid_pre))
# 
# grid_pts <- cbind(grid_pts, rep(psi_seq, times = nas))
# 
# colnames(grid_pts) <- c(paste0("alpha", 2:M), "psi")



# version that randomly generates alphas every time

num_psi <- length(psi_seq)

grid_pts <- matrix(NA, nrow = nas * num_psi, ncol = M)

ctr <- 1

for (sn in 1:nas) {
  
  for (sp in 1:num_psi) {
    
    grid_pts[ctr, 1:(M - 1)] <- rtruncnorm(n=M-1, a=alpha_range[1], b=alpha_range[2], mean=0, sd=alpha_sd) # runif(M - 1, min = alpha_range[1], max = alpha_range[2])
    
    grid_pts[ctr, M] <- psi_seq[sp]
    
    ctr <- ctr + 1
    
  }
  
}

colnames(grid_pts) <- c(paste0("alpha", 2:M), "psi")

sink_filename <- here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "ordered_grid_canon.txt")) 

########## Start of parallel loops to calculate canonical statistics
# Time required will exponentially increase with M.

canon_ests <- foreach(pt = iter(grid_pts, by = "row"), .combine = rbind, .packages = c("potts", "here"), .inorder = F, .options.RNG = seed) %dorng% {
  
  .GlobalEnv$x_dim <- x_dim
  .GlobalEnv$y_dim <- y_dim
  .GlobalEnv$M <- M
  .GlobalEnv$num_realizations <- num_realizations
  .GlobalEnv$record_progress <- record_progress
  .GlobalEnv$sink_filename <- sink_filename
  
  
  
  ordered_pt <- c(sort(c(0, pt[-M]), na.last = T), pt[M])
  
  sm_ests <- rep(NA, length = M)
  
  # Generate canonical statistics given the tuple of psi and alpha values given in grid_pts
  
  # Parameters of the Potts model
  theta <- ordered_pt
  x <- matrix(1, nrow = y_dim, ncol = x_dim)
  foo <- packPotts(x, M)
  
  out <- potts(foo, theta, nbatch = 200*num_realizations, blen = 1, boundary = "free")
  
  canon_stat <- out$batch[seq(200, 200*num_realizations, by = 200), ]
  
  sm_ests <- apply(canon_stat, 2, mean)
  
  if (record_progress) {
    sink(sink_filename, append=T)
    cat(paste0(paste(c(ordered_pt, sm_ests), collapse = " "), "\n"))
    sink()
  }
  
  sm_ests
  
}


