library(ggplot2)



# find info matrices grid_nbr_idx and grid_nbr_type
grid_info_mats <- function(num_gp = 130, num_spp, 
                           neighbors = getNeighbors(matrix(1, 13, 10), c(2, 2, 2, 2)), 
                           cluster_labs) {
  
  grid_nbr_idx <- array(NA, dim = c(num_gp, 9), 
                        dimnames = list(paste0("grid", 1:num_gp), 
                                        paste0("nbr", 0:8)))
  
  # convention: 
  # 0 is self
  # 1 is west neighbor
  # 2 is northwest neighbor
  # 3 is north neighbor
  # 4 is northeast neighbor
  # 5 is east
  # 6 is southeast
  # 7 is south
  # 8 is southwest
  
  # correspondence between convention above and the order of neighbor indices from getNeighbors
  nbr_order <- c(1, 8, 4, 6, 2, 7, 3, 5)
  
  
  
  for (j in 1:num_gp) {
    for (k in 1:9) {
      
      if (k == 1) {
        
        grid_nbr_idx[j, k] <- j
        
      } else {
        
        nbr_idx <- neighbors[j, nbr_order[k - 1]]
        
        if (nbr_idx != (num_gp + 1)) {
          grid_nbr_idx[j, k] <- nbr_idx
        }
        
      }
      
    }
  }
  
  
  
  hv_grid_idx <- c(1, 3, 5, 7)
  d_grid_idx <- c(2, 4, 6, 8)
  
  hv_idx <- which(c(0:8) %in% hv_grid_idx)
  d_idx <- which(c(0:8) %in% d_grid_idx)
  
  gp_nbr_type <- matrix(NA, nrow = num_spp, ncol = num_gp)
  
  for (i in 1:num_spp) {
    
    for (j in 1:num_gp) {
      
      nbr_labs <- cluster_labs[i, grid_nbr_idx[j, ]]
      
      # skips over any grid regions at the corner or on the edge
      if (any(is.na(nbr_labs))) {
        gp_nbr_type[i, j] <- "edge"
        next
      }
      
      if (nbr_labs[1] == 2) { # first considering the case that the center grid region is Thomas
        
        # notation: csr#hv#d, for number of horizontal/vertical and diagonal neighbors
        hv <- sum(which(nbr_labs == 1) %in% hv_idx)
        d <- sum(which(nbr_labs == 1) %in% d_idx)
        
        gp_nbr_type[i, j] <- paste0("csr", hv, "hv", d, "d")
        
      } else { # this is the most common case, i.e. all regions can be considered CSR
        gp_nbr_type[i, j] <- "csr9"
      }
      
    }
    
  }
  
  
  
  return(list(gp_nbr_type = gp_nbr_type, grid_nbr_idx = grid_nbr_idx))
  
}



# Set covariance function, given which grid region 0-8 I'm focusing on
# eps: single number indicating the precision in the form of pixel width/height
# subsetted to the angles with positive support
# precision and subsetting is for rotmean function
grid_region_setcov <- function(l = 0, eps = 0.01) {
  
  g0 <- matrix(0, nrow = 3, ncol = 3)
  g0[2, 2] <- 1
  
  gl <- matrix(0, nrow = 3, ncol = 3)
  if (l == 0) {
    gl[2, 2] <- 1
  } else if (l == 1) {
    gl[2, 1] <- 1
  } else if (l == 2) {
    gl[3, 1] <- 1
  } else if (l == 3) {
    gl[3, 2] <- 1
  } else if (l == 4) {
    gl[3, 3] <- 1
  } else if (l == 5) {
    gl[2, 3] <- 1
  } else if (l == 6) {
    gl[1, 3] <- 1
  } else if (l == 7) {
    gl[1, 2] <- 1
  } else if (l == 8) {
    gl[1, 1] <- 1
  }
  
  w <- as.im(owin(mask = g0 == 1), eps = eps)
  
  wl <- as.im(owin(mask = gl == 1), eps = eps)
  cw <- setcov(wl, w)
  
  # subsetting the set covariance accordingly
  if (l == 0) {
    subset_win <- owin(c(-1, 1), c(-1, 1))
  } else if (l == 1) {
    subset_win <- owin(c(-2, 0), c(-1, 1))
  } else if (l == 2) {
    subset_win <- owin(c(-2, 0), c(0, 2))
  } else if (l == 3) {
    subset_win <- owin(c(-1, 1), c(0, 2))
  } else if (l == 4) {
    subset_win <- owin(c(0, 2), c(0, 2))
  } else if (l == 5) {
    subset_win <- owin(c(0, 2), c(-1, 1))
  } else if (l == 6) {
    subset_win <- owin(c(0, 2), c(-2, 0))
  } else if (l == 7) {
    subset_win <- owin(c(-1, 1), c(-2, 0))
  } else if (l == 8) {
    subset_win <- owin(c(-2, 0), c(-2, 0))
  }
  
  cw_subset <- cw[subset_win]
  
  return(cw_subset)
  
}



# Function to return the isotropic translation edge correction function
# I arbitrarily set adjust = 200 to get a decent-looking rotational mean curve. 
# May need to check later
rigid_motion_correction <- function(cw_subset, adjust = 200) {
  
  # Taking advantage of padzero = FALSE to compute rotational mean only over the parts where the set covariance is nonzero
  # Origin is automatically at (0, 0), due to how grid_region_setcov is constructed
  CW <- rotmean(cw_subset, padzero = FALSE, adjust = adjust)
  
  return(as.function(CW))
  
}

# Corresponding functions for CSR scenario
# 
# The CSR pcf is just one. Simplifies quite a bit of things


csrpcf0 <- function(r) {rmc0(r)}
csrpcf1 <- function(r) {rmc1(r) * 0.5}
csrpcf2 <- function(r) {rmc2(r) * 0.25}


csrpcf0_times2pir <- function(r) {
  rmc0(r) * 2 * pi * r
}

csrK0 <- Vectorize(function(r) {
  
  integrate(csrpcf0_times2pir, lower = 0, upper = r)$value
  
}, vectorize.args = "r")

csrpcf1_times2pir <- function(r) {
  rmc1(r) * pi * r
}

csrK1 <- Vectorize(function(r) {
  
  integrate(csrpcf1_times2pir, lower = 0, upper = r)$value
  
}, vectorize.args = "r")

csrpcf2_times2pir <- function(r) {
  rmc2(r) * pi/2 * r
}

csrK2 <- Vectorize(function(r) {
  
  integrate(csrpcf2_times2pir, lower = 0, upper = r)$value
  
}, vectorize.args = "r")



# pcf for Thomas point process, Matern I point process, and LGCP

# r: evaluate pcf at this value
# kappa: intensity of the Poisson process of cluster centres.
# scale: Standard deviation of random displacement (along each coordinate axis) of a point from its cluster center
thomaspcf <- function(r, kappa, scale) {
  1 + exp(-r^2 / (4 * scale^2)) / (4 * pi * kappa * scale^2)
}

# r: evaluate pcf at this value
# kappa: intensity of the Poisson process of proposal points. A single positive number.
# delta: inhibition distance. Proposal points within distance delta from another proposal point are deleted.
maternIpcf <- Vectorize(function(r, kappa, delta) {
  
  if (r <= delta) { # technically, pcf is undefined at r = delta.
    return(0)
  } else {
    (kappa / exp(-kappa * pi * delta^2)) * exp(-kappa * circle_intersection(r, delta))
  }
  
}, vectorize.args = "r")

# # plot
# plot(r_grid, maternIpcf(r_grid, kappa = 2, delta = 0.4), type = "l")

exp_covfn <- function(r, var, scale) {
  var * exp(-(r / scale))
}

lgcppcf <- function(r, covfn = exp_covfn, var, scale) {
  exp(covfn(r, var, scale))
}



# calculate area of intersection of two circles with same radius delta, with centers separated by a distance r
# adapted from: https://stackoverflow.com/questions/44437793/how-to-calculate-the-intersection-area-of-two-circles-in-shiny-or-r-code
circle_intersection <- Vectorize(function(r, delta){
  delta_sq <- delta * delta
  
  if (r > delta + delta) # Circles do not overlap
  {
    return(0)
  } else { # Circles partially overlap
    phi <- (acos((delta_sq + (r * r) - delta_sq) / (2 * delta * r))) * 2
    theta <- (acos((delta_sq + (r * r) - delta_sq) / (2 * delta * r))) * 2
    area2 <- 0.5 * theta * delta_sq - 0.5 * delta_sq * sin(theta)
    area1 <- 0.5 * phi * delta_sq - 0.5 * delta_sq * sin(phi)
    return(area1 + area2)
  }
}, vectorize.args = "r")

# # testing
# circle_intersection(r = 2, delta = 2) # 4.913479
# 
# r <- 2
# r^2 * (2 * pi / 3 - sqrt(3) / 2) # alt formula for specific case that circles have same radius, and are radius apart from each other
# 
# # plot
# plot(r_grid, circle_intersection(r_grid, delta = 0.4), type = "l")



# mapping grid region index l = 0, ..., 8 to corresponding matrix indices

l2matidx <- Vectorize(function(l = 0) {
  
  mat_idx <- c(NA, NA)
  
  if (l == 0) {
    mat_idx[1] <- 2
    mat_idx[2] <- 2
  } else if (l == 1) {
    mat_idx[1] <- 2
    mat_idx[2] <- 1
  } else if (l == 2) {
    mat_idx[1] <- 3
    mat_idx[2] <- 1
  } else if (l == 3) {
    mat_idx[1] <- 3
    mat_idx[2] <- 2
  } else if (l == 4) {
    mat_idx[1] <- 3
    mat_idx[2] <- 3
  } else if (l == 5) {
    mat_idx[1] <- 2
    mat_idx[2] <- 3
  } else if (l == 6) {
    mat_idx[1] <- 1
    mat_idx[2] <- 3
  } else if (l == 7) {
    mat_idx[1] <- 1
    mat_idx[2] <- 2
  } else if (l == 8) {
    mat_idx[1] <- 1
    mat_idx[2] <- 1
  }
  
  return(mat_idx)
  
})



# Returns a corresponding column-wise index vector
mat_idx2vec <- function(row_idx, col_idx, y_dim) {
  
  if (any(row_idx > y_dim)) {stop("row_idx greater than y_dim")}
  
  return(row_idx + (col_idx - 1) * y_dim)
  
}



# assuming the vector is col-wise
vec_idx2mat <- function(vec_idx, y_dim) {
  
  return(list(row_idx = (vec_idx - 1) %% y_dim + 1, col_idx = (vec_idx - 1) %/% y_dim + 1))
  
}
# will be used to get number of cells in each grid region in each subject



# Return a window that's the union of specified grid regions (given a logical matrix).
# Assumes the (1,1) entry of the logical matrix corresponds to owin(xrange=c(0,1), yrange=c(0,1)),
# and the (i,j) entry of matrix corresponds to the ith row above and
# the jth column to the right of this origin square
partial_window <- function(mask) {
  
  specified_grids <- which(mask, arr.ind = T)
  
  grid_list <- vector("list", length = nrow(specified_grids))
  
  for (i in 1:nrow(specified_grids)) {
    
    mat_idx <- specified_grids[i, ]
    grid_list[[i]] <- owin(xrange = c(mat_idx[2] - 1, mat_idx[2]),
                           yrange = c(mat_idx[1] - 1, mat_idx[1]))
    
  }
  
  p_win <- do.call("union.owin", grid_list)
  
  return(p_win)
  
}



# Plotting the matrix. In contrast to the function plot.matrix, this function
# puts the [1, 1] entry at the lower left rather than the top left.
# This is to make plotting the cluster labels consistent with plotting the spatial point pattern.
plot_matrix = function(in_mat, main = "", discrete = T, levels = NULL){
  x <- seq(1, ncol(in_mat)) - 0.5
  y <- seq(1, nrow(in_mat)) - 0.5
  
  df <- expand.grid(X=x, Y=y)
  
  if (discrete) {
    if (is.null(levels)) {
      df$Z <- as.character(t(in_mat))
    } else {
      df$Z <- factor(t(in_mat), levels = levels)
    }
  } else {
    df$Z <- c(t(in_mat))
  }
  
  if (discrete) {
    ggplot(df, aes(X, Y, fill= Z)) + 
      geom_tile() +
      ggtitle(main) +
      scale_fill_viridis_d(name = "Cluster", drop = F) + 
      scale_x_continuous(breaks = x + 0.5) + 
      scale_y_continuous(breaks = y + 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  } else {
    ggplot(df, aes(X, Y, fill= Z)) + 
      geom_tile() +
      ggtitle(main) + 
      scale_fill_viridis_c(name = "Score") + 
      scale_x_continuous(breaks = x + 0.5) + 
      scale_y_continuous(breaks = y + 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
}



# Function to plot cluster labels with spatial point pattern superimposed
# Plotting the matrix. In contrast to the function plot.matrix, this function
# puts the [1, 1] entry at the lower left rather than the top left.
# This is to make plotting the cluster labels consistent with plotting the spatial point pattern.
plot_clust_spp = function(clust_spp, main = "", tile_x_incr = 0, tile_y_incr = 0, truncate = F, levels = NULL, plot_tile = T){
  in_mat <- clust_spp$cluster_lab
  
  X <- seq(1, ncol(in_mat)) - 0.5 + tile_x_incr
  Y <- seq(1, nrow(in_mat)) - 0.5 + tile_y_incr
  
  X_range <- c(min(X) - 0.5, max(X) + 0.5)
  Y_range <- c(min(Y) - 0.5, max(Y) + 0.5)
  
  df <- expand.grid(X=X, Y=Y)
  df$Z <- as.character(t(in_mat))
  
  spp <- clust_spp$spp
  df2 <- data.frame(x = spp$x, y = spp$y, Type = spp$marks)
  if (truncate) {
    df2 <- df2[df2$x > X_range[1] & df2$x < tail(X_range, 1) &
                 df2$y > Y_range[1] & df2$y < tail(Y_range, 1), ]  
  }
  
  if (is.null(levels)) {
    df$Z <- as.character(t(in_mat))
  } else {
    df$Z <- factor(t(in_mat), levels = levels)
  }
  
  if (plot_tile) {
    ggplot(df) + 
      geom_tile(aes(X, Y, fill= Z)) +
      scale_fill_viridis_d(name = "Cluster", drop = F) + 
      ggtitle(main) +
      scale_x_continuous(breaks = seq(1, ncol(in_mat))) + 
      scale_y_continuous(breaks = seq(1, nrow(in_mat))) + 
      geom_point(data = df2, mapping = aes(x, y, col = Type), size = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("")
  } else {
    ggplot(df) + 
      scale_fill_viridis_d(name = "Cluster", drop = F) + 
      ggtitle(main) +
      scale_x_continuous(breaks = seq(1, ncol(in_mat))) + 
      scale_y_continuous(breaks = seq(1, nrow(in_mat))) + 
      geom_point(data = df2, mapping = aes(x, y, col = Type), size = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("") + ylab("")
  }
  
}



# Function to generate cluster labels and associated spatial point pattern
generate_clusters_and_spp <- function(M, psi, S_dim = list(x_dim = 0:12, y_dim = 0:10),
                                      lambda_dot_vec = c(32, 20), prob_mat = cbind(c(0.5, 0.5), c(0.8, 0.2)), specs_idx = c(2, 2, 2, 1, 1)) {
  
  # TRUE BREAKPOINTS OF GRID REGIONS TO WHICH CLUSTERS WILL BE ASSIGNED
  x_dim_true <- S_dim$x_dim
  y_dim_true <- S_dim$y_dim
  
  num_gp <- (length(x_dim_true) - 1) * (length(y_dim_true) - 1)
  
  
  
  theta <- c(rep(0, M), psi)
  x <- matrix(1, nrow = length(y_dim_true) - 1, ncol = length(x_dim_true) - 1)
  foo <- packPotts(x, M)
  
  
  
  out <- potts(foo, theta, nbatch = 10, blen = 20, boundary = "free")
  
  cluster_lab <- unpackPotts(out$final)
  
  
  
  # Defining the pcfs
  # Thomas parameters
  thomas_kappa <- 5
  thomas_scale <- 0.1
  
  # Matérn I parameters (depends on lambda)
  intens_root <- function(kappa, delta, intens) {kappa * exp(-kappa * pi * delta^2) - intens}
  
  delta1 <- 0.037
  maternI_kappa1 <- uniroot(intens_root, interval = c(0, 1), delta1, intens = lambda_dot_vec[specs_idx[3]], extendInt = "upX")$root
  
  if (M >= 4) {
    # Second set of Thomas parameters 
    thomas_kappa2 <- 6
    thomas_scale2 <- 0.055
    
    if (M == 5) {
      # Matérn I parameters (depends on lambda)
      delta2 <- 0.02
      maternI_kappa2 <- uniroot(intens_root, interval = c(0, 1), delta2, intens = lambda_dot_vec[specs_idx[5]], extendInt = "upX")$root
    }
  }
  
  
  
  if (1 %in% c(cluster_lab)) {
    # CSR points
    w <- partial_window(mask = (cluster_lab == 1))
    csr_ppp <- rpoispp(lambda = lambda_dot_vec[specs_idx[1]], win = w)
    mark <- sample(c(1, 2), size = npoints(csr_ppp), replace = T, prob = prob_mat[, specs_idx[1]])
  } else {
    csr_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(csr_ppp) <- factor(mark, levels = c(1, 2))
  
  if (2 %in% c(cluster_lab)) {
    # Thomas process points
    w <- partial_window(mask = (cluster_lab == 2))
    thomas_ppp <- rThomas(kappa = thomas_kappa, scale = thomas_scale, mu = lambda_dot_vec[specs_idx[2]] / thomas_kappa, win = w)
    mark <- sample(c(1, 2), size = npoints(thomas_ppp), replace = T, prob = prob_mat[, specs_idx[2]])
  } else{
    thomas_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(thomas_ppp) <- factor(mark, levels = c(1, 2))
  
  if (3 %in% c(cluster_lab)) {
    # Matern I process points
    w <- partial_window(mask = (cluster_lab == 3))
    maternI_ppp <- rMaternI(kappa = maternI_kappa1, r = delta1, win = w)
    mark <- sample(c(1, 2), size = npoints(maternI_ppp), replace = T, prob = prob_mat[, specs_idx[3]])
  } else{
    maternI_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(maternI_ppp) <- factor(mark, levels = c(1, 2))
  
  if (4 %in% c(cluster_lab)) {
    # Thomas process points, 2nd set
    w <- partial_window(mask = (cluster_lab == 4))
    thomas_ppp2 <- rThomas(kappa = thomas_kappa2, scale = thomas_scale2, mu = lambda_dot_vec[specs_idx[4]] / thomas_kappa2, win = w)
    mark <- sample(c(1, 2), size = npoints(thomas_ppp2), replace = T, prob = prob_mat[, specs_idx[4]])
  } else{
    thomas_ppp2 <- ppp() # empty point pattern
    mark <- c()
  }
  marks(thomas_ppp2) <- factor(mark, levels = c(1, 2))
  
  if (5 %in% c(cluster_lab)) {
    # Matern I process points
    w <- partial_window(mask = (cluster_lab == 5))
    maternI_ppp2 <- rMaternI(kappa = maternI_kappa2, r = delta2, win = w)
    mark <- sample(c(1, 2), size = npoints(maternI_ppp2), replace = T, prob = prob_mat[, specs_idx[5]])
  } else{
    maternI_ppp2 <- ppp() # empty point pattern
    mark <- c()
  }
  marks(maternI_ppp2) <- factor(mark, levels = c(1, 2))
  
  spp <- superimpose(csr_ppp, thomas_ppp, thomas_ppp2, maternI_ppp, maternI_ppp2, 
                     W = owin(xrange = c(x_dim_true[1], tail(x_dim_true, 1)), yrange = c(y_dim_true[1], tail(y_dim_true, 1))))
  
  
  
  return(list(cluster_lab = cluster_lab, spp = spp))
  
}





# Function to generate spp's pixel by pixel, instead of by cluster to ensure true conditional independence.
# For the M = 5 scenario.
generate_spp_by_pixel <- function(cluster_im, x_dim_true, y_dim_true, lambda, thomas_kappa, thomas_scale, lgcp_mu, lgcp_var, lgcp_scale, thomas_kappa2, thomas_scale2, maternI_kappa, delta) {
  
  num_gp <- length(c(cluster_im))
  
  pixel_spp_list <- vector(mode = "list", length = num_gp)
  
  mask_mat <- matrix(FALSE, nrow = nrow(cluster_im), ncol = ncol(cluster_im))
  
  for (k in 1:num_gp) {
    
    # focusing on kth pixel
    mask_mat[k] <- TRUE
    
    if (cluster_im[k] == 1) {
      # CSR points
      w <- partial_window(mask = mask_mat)
      pixel_spp_list[[k]] <- rpoispp(lambda = lambda, win = w)
    } else if (cluster_im[k] == 2) {
      w <- partial_window(mask = mask_mat)
      pixel_spp_list[[k]] <- rThomas(kappa = thomas_kappa, scale = thomas_scale, mu = lambda / thomas_kappa, win = w)
    } else if (cluster_im[k] == 3) {
      w <- partial_window(mask = mask_mat)
      pixel_spp_list[[k]] <- rLGCP(model = "exp", mu = lgcp_mu, var = lgcp_var, scale = lgcp_scale, win = w)
    } else if (cluster_im[k] == 4) {
      w <- partial_window(mask = mask_mat)
      pixel_spp_list[[k]] <- rThomas(kappa = thomas_kappa2, scale = thomas_scale2, mu = lambda / thomas_kappa2, win = w)
    } else if(cluster_im[k] == 5) {
      w <- partial_window(mask = mask_mat)
      pixel_spp_list[[k]] <- rMaternI(kappa = maternI_kappa, r = delta, win = w)
    }
    
    # resetting the mask
    mask_mat[k] <- FALSE
    
  }
  
  X <- do.call("superimpose", pixel_spp_list)
  
  # also superimposing X with empty window to ensure 
  # window(X) contains all the unit squares across the subject
  X <- superimpose(X, 
                   ppp(window = owin(xrange = c(x_dim_true[1], tail(x_dim_true, 1)), yrange = c(y_dim_true[1], tail(y_dim_true, 1)))))
  
  return(X)
  
}



## Whole pipeline from retrieving the canonical statistics and grid points to integration fitting
# @input parms: referring to parameters other than the first alpha offset
# @other_mods: assuming this is a vector of at least one element
# @other_mod_parms: 
potts_interp_testing <- function(M, x_dim_length, y_dim_length, parms,
                                 psi_seq = c(seq(0, 2.4, by = 0.2), 2.5), 
                                 alpha_seq = seq(-5, 5, by = 0.05), 
                                 psi_test = T,
                                 alpha_test = T,
                                 other_mods = c("xgboost"), 
                                 other_mod_parms = vector(mode = "list", length = length(other_mods))) {
  
  # Take 251 points from the knn, fit a 10th-degree smoothing spline, then calculate the quantities
  
  ptm_all <- proc.time()
  
  res_string <- paste0("res", y_dim_length, "x", x_dim_length)
  
  num_cores=detectCores()
  cl <- makeCluster(num_cores - 2) 
  registerDoParallel(cl)
  
  memory_vec <- rep(NA, length(other_mods))
  names(memory_vec) <- other_mods
  
  MSE_mat <- matrix(NA, nrow = 2, ncol = length(other_mods))
  colnames(MSE_mat) <- other_mods
  
  # Setting up the knnreg interpolation and other_mods
  
  if ("knnreg" %in% other_mods | "xgboost" %in% other_mods) {
    
    if ("knnreg" %in% other_mods) {
      # Reading from .rds file
      potts_samps <- readRDS(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "potts_samps.rds")))
      knn_mod_alpha_psi <- append(potts_samps[[2]], list(potts_samps[[1]]))
      
      obj_size <- object.size(knn_mod_alpha_psi)
      # print(paste0("knnreg: ", obj_size, " bytes"))
      memory_vec[which(other_mods == "knnreg")] <- obj_size
    }
    
    if ("xgboost" %in% other_mods) {
      
      if (isTRUE(other_mod_parms[[which(other_mods == "xgboost")]])) {
        # read in original canonical statistics
        ordered_combined <- as.data.frame(read_delim(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "ordered_grid_canon.txt")), delim = " ", show_col_types = FALSE))
        grid_pts_ordered <- ordered_combined[, 1:(M+1)]
        canon_vals_ordered <- ordered_combined[, (M+2):(2*(M + 1))]
        colnames(grid_pts_ordered) <- c(paste0("alpha", 1:M), "psi")
        
        xgboost_mod_list <- vector(mode = "list", length = M+1)
        param <- list(booster = "gbtree", eval_metric = "rmse")
        for (i in 1:(M+1)) {
          dtrain <- xgb.DMatrix(data = as.matrix(grid_pts_ordered), label = canon_vals_ordered[, i])
          xgboost_mod_list[[i]] <- xgb.train(param, dtrain, nrounds=100, max_depth = 12)
        }
      } else {
        # read in file made through potts_samps_xgboost.R instead
        potts_samps_xgb <- readRDS(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "potts_samps_xgb.rds")))
        xgboost_mod_list <- append(potts_samps_xgb[[2]], list(potts_samps_xgb[[1]]))
      }
      
      obj_size <- object.size(xgboost_mod_list)
      # print(paste0("xgboost: ", obj_size, " bytes"))
      memory_vec[which(other_mods == "xgboost")] <- obj_size
      
    }
    
  }
  
  if ("randomforest" %in% other_mods | "mars" %in% other_mods) {
    
    # read in original canonical statistics
    ordered_combined <- as.data.frame(read_delim(here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "ordered_grid_canon.txt")), delim = " ", show_col_types = FALSE))
    grid_pts_ordered <- ordered_combined[, 1:(M+1)]
    canon_vals_ordered <- ordered_combined[, (M+2):(2*(M + 1))]
    colnames(grid_pts_ordered) <- c(paste0("alpha", 1:M), "psi")
    
    if ("randomforest" %in% other_mods) {
      # these model lists are ordered from the model for the minimum, model for the 2nd smallest, to the model for the maximum
      rf_mod_list <- vector(mode = "list", length = M+1)
      for (i in 1:(M+1)) {
        rf_mod_list[[i]] <- ranger(y = canon_vals_ordered[, i], x = grid_pts_ordered)
      }
      # print(paste0("randomforest: ", object.size(rf_mod_list), "bytes"))
    }
    
    if ("mars" %in% other_mods) {
      mars_mod_list <- vector(mode = "list", length = M+1)
      for (i in 1:(M+1)) {
        mars_mod_list[[i]] <- earth(y = canon_vals_ordered[, i], x = grid_pts_ordered, 
                                    degree = 3)
      }
      # print(paste0("mars: ", object.size(mars_mod_list), "bytes"))
    }
    
  }
  
  
  
  # psi sequence to integrate along
  psi_seq_inter <- c(seq(0, 2.4, by = 0.2), 2.5)
  psi_incr_vec <- psi_seq_inter[-length(psi_seq_inter)] 
  psi_incr <- psi_incr_vec[2] - psi_incr_vec[1] # assuming an equal increment
  psi_bins <- psi_incr_vec + psi_incr/2 # seq(0.1, 2.5, by = 0.2)
  
  # alpha sequence to integrate along
  alpha_seq_inter = seq(-5, 5, by = 0.2)
  alpha_incr <- (range(alpha_seq_inter)[2] - range(alpha_seq_inter)[1]) / 10 # assuming an equal increment
  alpha_incr_vec <- sort(c(seq(-5 + alpha_incr/2, 5 - alpha_incr/2, by = alpha_incr), 0))
  alpha_bins <- alpha_incr_vec + alpha_incr/2
  alpha_bins[alpha_bins == 0] <- -0.5 
  
  
  
  # # above psi_seq_inter and alpha_seq_inter are just to get the psi_bins and alpha_bins. Here are the actual sequences for interpolation:
  # psi_seq <- c(seq(0, 2.4, by = 0.2), 2.5)
  # alpha_seq <- seq(-5, 5, by = 0.05)
  # # commented out because I use function inputs instead
  
  
  
  if (psi_test) {
    
    if ("knnreg" %in% other_mods) {
      buffer_idx <- 3
      mod_idx <- seq(1, buffer_idx + length(other_mods))
      mod_vec <- c("benchmark", "raw knnreg", "smoothed knnreg", other_mods)
    } else {
      buffer_idx <- 1
      mod_idx <- seq(1, buffer_idx + length(other_mods))
      mod_vec <- c("benchmark", other_mods)
    }
    
    # benchmark estimates for the given parms (last element dropped, replaced with psi_seq)
    psi_sequence <- seq(0, 2.5, by = 0.05)
    psi_sequence_mat <- cbind(rep(0, length(psi_sequence)), matrix(parms[-M], nrow = length(psi_sequence), ncol = M - 1, byrow = T), psi_sequence)
    benchmark_ests_psi <- foreach(i = 1:nrow(psi_sequence_mat), .combine = rbind, .packages = c("potts", "here"), .options.RNG = 942) %dorng% {
      
      num_realizations <- 1000
      
      
      
      pt <- psi_sequence_mat[i, ]
      
      sm_ests <- rep(NA, length = M)
      
      # Generate canonical statistics given the tuple of psi and alpha values given in psi_sequence_mat
      
      # Parameters of the Potts model
      theta <- pt
      x <- matrix(1, nrow = y_dim_length, ncol = x_dim_length)
      foo <- packPotts(x, M)
      
      out <- potts(foo, theta, nbatch = 200*num_realizations, blen = 1, boundary = "free")
      
      canon_stat <- out$batch[seq(200, 200*num_realizations, by = 200), ]
      
      sm_ests <- apply(canon_stat, 2, mean)
      
      sm_ests
      
    }
    plot(psi_sequence, benchmark_ests_psi[, ncol(benchmark_ests_psi)], type = "l", col = 1, xlab = expression(psi), ylab = "Expected Number of Matching Adjacent Pairs")
    
    psi_linear_interp <- approxfun(psi_sequence, benchmark_ests_psi[, 2])
    lin_interp <- psi_linear_interp(psi_seq)
    
    
    
    # for psi parameter
    
    # assembling matrix for prediction
    # in the order of other alphas (if any) and psi
    preds_psi <- rep(NA, length(psi_seq))
    
    psi_seq_mat <- cbind(rep(0, length(psi_seq)), matrix(parms[-M], nrow = length(psi_seq), ncol = M - 1, byrow = T), psi_seq)
    colnames(psi_seq_mat) <- c(paste0("alpha", 1:M), "psi")
    
    sm <- M + 1
    
    if ("randomforest" %in% other_mods) {
      preds_rf <- rep(NA, length(psi_seq))
    }
    
    if ("mars" %in% other_mods) {
      preds_mars <- rep(NA, length(psi_seq))
    }
    
    if ("xgboost" %in% other_mods) {
      preds_xgboost <- rep(NA, length(psi_seq))
    }
    
    for (i in 1:length(psi_seq)) {
      pt <- as.numeric(psi_seq_mat[i, ])
      offsets_order <- c(order(pt[1:M]), M+1)
      pt_ordered <- psi_seq_mat[i, offsets_order, drop = F]
      colnames(pt_ordered) <- c(paste0("alpha", 1:M), "psi")
      
      if ("knnreg" %in% other_mods) {
        incr_vec <- rep(NA, length = M+1)
        for (incr_sm in 1:(M+1)) {
          if (incr_sm <= M) {
            incr_vec[incr_sm] <- which(pt_ordered[incr_sm] <= alpha_bins)[1]
          } else {
            incr_vec[incr_sm] <- which(pt_ordered[incr_sm] <= psi_bins)[1]
          }
        }
        preds_psi[i] <- predict(knn_mod_alpha_psi[[M+1]][[incr_vec]], pt_ordered)  
      }
      
      # Other modeling
      if ("randomforest" %in% other_mods) {
        preds_rf[i] <- predict(rf_mod_list[[M+1]], data.frame(pt_ordered))$predictions
      }
      
      if ("mars" %in% other_mods) {
        preds_mars[i] <- as.numeric(predict(mars_mod_list[[M+1]], data.frame(pt_ordered)))
      }
      
      if ("xgboost" %in% other_mods) {
        xgb_test = xgb.DMatrix(data = pt_ordered)
        preds_xgboost[i] <- predict(xgboost_mod_list[[M+1]], xgb_test)
      }
      
    }
    
    
    
    if ("knnreg" %in% other_mods) {
      
      lines(psi_seq, preds_psi, type = "l", col = 2)
      
      # for integration for psi
      
      p <- 10
      
      # for psi spline calculations for log(d(theta))
      PSI <- matrix(psi_seq, length(psi_seq), p, byrow = F) 
      pow <- matrix(1:p, length(psi_seq), p, byrow = T)
      X <- PSI^pow
      beta <- lm(preds_psi ~ X)$coef # replace psi_samps with the predictions
      
      lines(psi_seq, predict(lm(preds_psi ~ X)), col = 3)
      
      MSE <- mean(preds_psi - lin_interp)^2
      # print(paste0("knnreg psi MSE: ", MSE))
      MSE_mat[1, which(other_mods == "knnreg")] <- MSE
      
    }
    
    # Other modeling
    if ("randomforest" %in% other_mods) {
      lines(psi_seq, preds_rf, col = buffer_idx + which(other_mods == "randomforest"), type = "l")
      # print(paste0("randomforest psi MSE: ", mean(preds_rf - lin_interp)^2))
    }
    
    if ("mars" %in% other_mods) {
      lines(psi_seq, preds_mars, col = buffer_idx + which(other_mods == "mars"), type = "l")
      # print(paste0("mars psi MSE: ", mean(preds_mars - lin_interp)^2))
    }
    
    if ("xgboost" %in% other_mods) {
      lines(psi_seq, preds_xgboost, col = buffer_idx + which(other_mods == "xgboost"), type = "l")
      
      MSE <- mean(preds_xgboost - lin_interp)^2
      # print(paste0("xgboost psi MSE: ", MSE))
      MSE_mat[1, which(other_mods == "xgboost")] <- MSE
    }
    
    legend("bottomright", lty = 1, col = mod_idx, legend = mod_vec)
    
    
    
    psi_integ <- Vectorize(function(psi_prime) {integrate(approxfun(psi_sequence, benchmark_ests_psi[, ncol(benchmark_ests_psi)]), range(psi_sequence)[1], psi_prime, rel.tol=0.1)$value})
    
    preds_psi_integ <- psi_integ(psi_sequence) # keep in mind that I'll just need to calculate 2 values at a time, instead of a sequence like this
    
    # plot(psi_sequence, preds_psi_integ, type = "l", main = "Differences ", ylab= "Integration values")
    # 
    # if ("knnreg" %in% other_mods) {
    #   
    #   pow <- matrix(1:p, length(psi_seq), p, byrow = T)
    #   
    #   X_integ <- cbind(psi_seq, psi_seq^(pow + 1) / (pow + 1))
    #   
    #   lines(psi_seq, X_integ %*% beta, type = "l", col = 2)
    #   
    #   legend("topleft", lty = 1, col = 1:2, legend = c("Benchmark Integration", 
    #                                                    "Polynomial Integration"))
    #   
    # }
    # 
    # if ("xgboost" %in% other_mods) {
    #   psi_integ2 <- Vectorize(function(psi_prime) {integrate(approxfun(psi_seq, preds_xgboost), range(psi_seq)[1], psi_prime, rel.tol=0.1)$value})
    #   lines(psi_seq, psi_integ2(psi_seq), col = 2)
    #   
    #   legend("topleft", lty = 1, col = 1:2, legend = c("Benchmark Integration", 
    #                                                    "XGBoost Integration"))
    # }
    
  }
  
  
  
  if (alpha_test) {
    
    if ("knnreg" %in% other_mods) {
      buffer_idx <- 2
      mod_idx <- seq(1, buffer_idx + length(other_mods))
      mod_vec <- c("benchmark", "raw knnreg", other_mods)  
    } else {
      buffer_idx <- 1
      mod_idx <- seq(1, buffer_idx + length(other_mods))
      mod_vec <- c("benchmark", other_mods)  
    }
    
    # benchmark estimates for the above parms (first element dropped, replaced with alpha_sequence, alpha1 still zero)
    alpha_sequence <- seq(-5, 5, by = 0.05)
    alpha_sequence_mat <- cbind(rep(0, length(alpha_sequence)), alpha_sequence, matrix(parms[-1], nrow = length(alpha_sequence), ncol = M - 1, byrow = T))
    benchmark_ests_alpha <- foreach(i = 1:nrow(alpha_sequence_mat), .combine = rbind, .packages = c("potts", "here"), .options.RNG = 942) %dorng% {
      
      .GlobalEnv$M <- M
      .GlobalEnv$x_dim_length <- x_dim_length
      .GlobalEnv$y_dim_length <- y_dim_length
      .GlobalEnv$alpha_sequence_mat <- alpha_sequence_mat
      
      num_realizations <- 1000
      
      
      
      pt <- alpha_sequence_mat[i, ]
      
      sm_ests <- rep(NA, length = M)
      
      # Generate canonical statistics given the tuple of psi and alpha values given in alpha_sequence_mat
      
      # Parameters of the Potts model
      theta <- pt
      x <- matrix(1, nrow = y_dim_length, ncol = x_dim_length)
      foo <- packPotts(x, M)
      
      out <- potts(foo, theta, nbatch = 200*num_realizations, blen = 1, boundary = "free")
      
      canon_stat <- out$batch[seq(200, 200*num_realizations, by = 200), ]
      
      sm_ests <- apply(canon_stat, 2, mean)
      
      sm_ests
      
    }
    plot(alpha_sequence, benchmark_ests_alpha[, 2], type = "l", col = 1, xlab = expression(alpha), ylab = expression(paste("Expected Number of Cluster Labels ", eta)))
    
    alpha_linear_interp <- approxfun(alpha_sequence, benchmark_ests_alpha[, 2])
    lin_interp <- alpha_linear_interp(alpha_seq)
    
    # for the alpha parameters
    # make predictions for alphas
    
    preds_alpha <- rep(NA, length(alpha_seq))
    # combining alpha_seq with parms
    alpha_seq_mat <- cbind(rep(0, length(alpha_seq)), alpha_seq, matrix(parms[-1], nrow = length(alpha_seq), ncol = M - 1, byrow = T))
    colnames(alpha_seq_mat) <- c(paste0("alpha", 1:M), "psi")
    
    if ("randomforest" %in% other_mods) {
      preds_rf <- rep(NA, length(alpha_seq))
    }
    
    if ("mars" %in% other_mods) {
      preds_mars <- rep(NA, length(alpha_seq))
    }
    
    if ("xgboost" %in% other_mods) {
      preds_xgboost <- rep(NA, length(alpha_seq))
    }
    
    # using the recursive models
    
    sm <- 2
    
    for (i in 1:length(alpha_seq)) {
      pt <- as.numeric(alpha_seq_mat[i, ])
      offsets_order <- c(order(pt[1:M]), M+1)
      pt_ordered <- alpha_seq_mat[i, offsets_order, drop = F]
      colnames(pt_ordered) <- c(paste0("alpha", 1:M), "psi")
      
      if ("knnreg" %in% other_mods) {
        incr_vec <- rep(NA, length = M+1)
        for (incr_sm in 1:(M+1)) {
          if (incr_sm <= M) {
            incr_vec[incr_sm] <- which(pt_ordered[incr_sm] <= alpha_bins)[1]
          } else {
            incr_vec[incr_sm] <- which(pt_ordered[incr_sm] <= psi_bins)[1]
          }
        }
        
        preds_alpha[i] <- predict(knn_mod_alpha_psi[[which(offsets_order == sm)]][[incr_vec]], pt_ordered)
        if (i > 1) {
          if (preds_alpha[i] < preds_alpha[i - 1]) {
            preds_alpha[i] <- preds_alpha[i - 1]
          }
        }
      }
      
      if ("randomforest" %in% other_mods) {
        preds_rf[i] <- predict(rf_mod_list[[which(offsets_order == sm)]], data.frame(pt_ordered))$predictions
      }
      
      if ("mars" %in% other_mods) {
        preds_mars[i] <- as.numeric(predict(mars_mod_list[[which(offsets_order == sm)]], data.frame(pt_ordered)))
      }
      
      if ("xgboost" %in% other_mods) {
        xgb_test = xgb.DMatrix(data = pt_ordered)
        
        preds_xgboost[i] <- predict(xgboost_mod_list[[which(offsets_order == sm)]], xgb_test)
        if (i > 1) {
          if (preds_xgboost[i] < preds_xgboost[i - 1]) {
            preds_xgboost[i] <- preds_xgboost[i - 1]
          }
        }
      }
    }
    
    if ("knnreg" %in% other_mods) {
      lines(alpha_seq, preds_alpha, type = "l", col = 2)
      
      MSE <- mean(preds_alpha - lin_interp)^2
      # print(paste0("knnreg alpha MSE: ", MSE))
      MSE_mat[2, which(other_mods == "knnreg")] <- MSE
    }
    
    if ("randomforest" %in% other_mods) {
      lines(alpha_seq, preds_rf, col = buffer_idx + which(other_mods == "randomforest"), type = "l")
      # print(paste0("randomforest alpha MSE: ", mean(preds_rf - lin_interp)^2))
    }
    
    if ("mars" %in% other_mods) {
      lines(alpha_seq, preds_mars, col = buffer_idx + which(other_mods == "mars"), type = "l")
      # print(paste0("mars alpha MSE: ", mean(preds_mars - lin_interp)^2))
    }
    
    if ("xgboost" %in% other_mods) {
      lines(alpha_seq, preds_xgboost, col = buffer_idx + which(other_mods == "xgboost"), type = "l")
      
      MSE <- mean(preds_xgboost - lin_interp)^2
      # print(paste0("xgboost alpha MSE: ", MSE))
      MSE_mat[2, which(other_mods == "xgboost")] <- MSE
    }
    
    legend("bottomright", lty = 1, col = mod_idx, legend = mod_vec)
    
    
    
    alpha_seq <- seq(-5, 5, by = 0.05)
    alpha_integ <- Vectorize(function(alpha_prime) {integrate(approxfun(alpha_sequence, benchmark_ests_alpha[, 2]), range(alpha_sequence)[1], alpha_prime, rel.tol=0.1)$value})
    
    preds_alpha_integ <- alpha_integ(alpha_sequence) # keep in mind that I'll just need to calculate 2 values at a time, instead of a sequence like this
    
    # plot(alpha_sequence, preds_alpha_integ, type = "l", main = "Differences ", ylab= "Integration values")
    # 
    # if ("xgboost" %in% other_mods) {
    #   alpha_integ2 <- Vectorize(function(alpha_prime) {integrate(approxfun(alpha_seq, preds_xgboost), range(alpha_seq)[1], alpha_prime, rel.tol=0.1)$value})
    #   lines(alpha_seq, alpha_integ2(alpha_seq), col = 2)
    #   
    #   legend("topleft", lty = 1, col = 1:2, legend = c("Benchmark Integration", 
    #                                                    "XGBoost Integration"))
    # }
    
  }
  
  # # deleting memory-heavy object(s). Not sure if this is necessary.
  # if ("randomforest" %in% other_mods) {
  #   rm(rf_mod_list)
  # }
  # 
  # if ("mars" %in% other_mods) {
  #   rm(mars_mod_list)
  # }
  # 
  # if ("xgboost" %in% other_mods) {
  #   rm(xgboost_mod_list)
  # }
  
  print(proc.time() - ptm_all)
  
  stopCluster(cl)
  
  return(list(memory_vec = memory_vec, MSE_mat = MSE_mat))
  
}


