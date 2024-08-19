# plotting helper functions

# Function to return data frame of cell locations and cell type, data for a single subject in wide format
# assumes types indicates cells of unrecorded type as "Other"
phen_df_format <- function(dat_i, types) {
  
  phen_names <- types[types != "Other"]
  
  # make a list of dataframes, one for each phenotype
  phen_list <- vector("list", length = length(phen_names) + 1)
  
  # cell centroid coordinates
  s <- dat_i[, 1:2]
  
  for (phen in 1:(length(phen_names) + 1)) {
    
    if (phen == (length(phen_names) + 1)) { # Indicating rows with no cell type, i.e. missing type
      type_missing <- rowSums(dat_i[, -c(1:2)]) == 0
      phen_list[[(length(phen_names) + 1)]] <- data.frame(s[type_missing,], Type = rep("Other", sum(type_missing)))
      next
    }
    
    phen_indicator <- dat_i[, colnames(dat_i) == phen_names[phen]]==1
    
    if (sum(phen_indicator) > 0) {
      phen_list[[phen]] <- data.frame(s[phen_indicator, ], Type = phen_names[phen])
    }
    
  }
  
  phen_pt_df <- do.call("rbind", args = c(phen_list, make.row.names = F))
  
  phen_pt_df$Type <- factor(phen_pt_df$Type, levels = types)
  
  return(phen_pt_df)
  
}



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Function to plot the spatial point pattern in data frame format
plot_df_spp <- function(phen_pt_df, main = "", legendTF = T, plotOther = T, 
                        type_labels = levels(phen_pt_df$Type), alpha = 0.3) {
  
  num_types <- length(levels(phen_pt_df$Type))
  col_vec <- gg_color_hue(num_types - 1)
  
  if (!plotOther) {
    phen_pt_df <- phen_pt_df[phen_pt_df$Type == "Other", ]
  } else {
    col_vec <- c(col_vec, "gray")
  }
  
  # Showing all cell phenotypes
  p <- ggplot(data = phen_pt_df, aes(x = Cell.X.Position, y = Cell.Y.Position, color = Type)) + 
    geom_point(alpha = alpha, show.legend = legendTF) +
    scale_color_manual(labels = type_labels, values = col_vec, drop = FALSE) +
    xlab(NULL) + ylab(NULL) +
    guides(color = guide_legend(override.aes = list(alpha = 1) ) ) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  if (main != "") {
    p <- p + ggtitle(main)
  }
  
  return(p)
  
}



# Function to plot cluster labels with real data spatial point pattern superimposed
# some differences from the original plot_clust_spp(): takes in the S_dim parameter instead of assuming that it's something like x_dim = 0:13, y_dim = 0:10. 
# Would automatically adjust for non-unit square tiles
plot_clust_spp_rd = function(clust_spp, S_dim, main = "", tile_x_incr = 0, tile_y_incr = 0, truncate = F, levels = NULL, types, type_labels = NULL){
  
  in_mat <- clust_spp$cluster_lab
  
  # getting the width and height of the tiles
  tile_width <- S_dim$x_dim[2] - S_dim$x_dim[1]
  tile_height <- S_dim$y_dim[2] - S_dim$y_dim[1]
  
  # getting the centers
  X <- S_dim$x_dim[-1] - tile_width/2 + tile_x_incr
  Y <- S_dim$y_dim[-1] - tile_height/2 + tile_y_incr
  
  X_range <- c(S_dim$x_dim[1], S_dim$x_dim[length(S_dim$x_dim)])
  Y_range <- c(S_dim$y_dim[1], S_dim$y_dim[length(S_dim$y_dim)])
  
  df <- expand.grid(X=X, Y=Y)
  
  if (is.null(levels)) {
    df$Z <- as.character(t(in_mat))
  } else {
    df$Z <- factor(t(in_mat), levels = levels)
  }
  
  spp <- clust_spp$spp
  df2 <- data.frame(x = spp$x, y = spp$y, Type = spp$marks)
  if (truncate) {
    df2 <- df2[df2$x > X_range[1] & df2$x < tail(X_range, 1) &
                 df2$y > Y_range[1] & df2$y < tail(Y_range, 1), ]  
  }
  
  if (!is.null(type_labels)) {
    df2$Type <- factor(df2$Type, levels = types)
    
    num_types <- length(type_labels)
    col_vec <- gg_color_hue(num_types - 1)
    
    col_vec <- c(col_vec, "gray")
  }
  
  ggplot(df) + 
    geom_tile(aes(X, Y, fill= Z)) +
    scale_fill_viridis_d(name = "Cluster", drop = F) + 
    ggtitle(main) +
    geom_point(data = df2, mapping = aes(x, y, col = Type), size = 0.5) + 
    scale_color_manual(labels = type_labels, values = col_vec, drop = FALSE) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
}



#' Iterate through the grid regions of a multitype spatial point pattern (starting from the one at the origin), finding 
#' the local intensities within them
#' @param S spatial point pattern of class ppp
#' @param S_dim list of grid region breakpoints on S's grid: list(x_dim = breakpoints of the grid regions on x-axis, y_dim = breakpoints of the grid regions on y-axis). 
#' Assumes breakpoints are equally spaced for each of x_dim and y_dim, so the grid regions are equally-sized rectangles
#' @param set0_as_missing if intensity for all types is 0 at a grid region, set intensities of all types to missing
#' @return list of local window pcfs, one for each grid region
local_window_calc_intens <- function(S, S_dim, set0_as_missing = F) {
  
  # iterate through the grid regions (starting from the one at the origin), finding 
  # the local pcf's within them
  
  ctr <- 1
  
  x_dim <- S_dim$x_dim
  y_dim <- S_dim$y_dim
  num_gp <- (length(x_dim) - 1) * (length(y_dim) - 1)
  
  mark <- marks(S)
  types <- levels(mark)
  if (!is.null(mark)) {
    local_intens <- matrix(NA, nrow = num_gp, ncol = length(types))
  } else{
    local_intens <- rep(NA, length = num_gp)
  }
  
  local_area <- spatstat.geom::area(owin(xrange = c(x_dim[1], x_dim[2]),
                                         yrange = c(y_dim[1], y_dim[2])))
  
  for (j in 1:(length(x_dim) - 1)) { 
    
    for (i in 1:(length(y_dim) - 1)) {
      
      local_domain <- owin(xrange = c(x_dim[j], x_dim[j + 1]),
                           yrange = c(y_dim[i], y_dim[i + 1]))
      
      # clipping S onto the window local_domain, omitting all points outside it
      S_subset <- S
      Window(S_subset) <- local_domain
      
      if (!is.null(mark)) {
        for (st in 1:length(types)) {
          local_intens[ctr, st] <- sum(marks(S_subset) == types[st], na.rm = T) / local_area
          
          # If local_domain is partially/fully outside window or there are no cells and set0_as_missing == T, set intensity as missing
          if ((set0_as_missing & npoints(S_subset) == 0) | !is.subset.owin(local_domain, Window(S))) {local_intens[ctr, st] <- NA}
        }
      } else {
        local_intens[ctr] <- npoints(S_subset) / local_area
        
        # If local_domain is partially/fully outside window or there are no cells and set0_as_missing == T, set intensity as missing
        if ((set0_as_missing & npoints(S_subset) == 0) | !is.subset.owin(local_domain, Window(S))) {local_intens[ctr] <- NA}
      }
      
      ctr <- ctr + 1
      
    }
    
  }
  
  return(local_intens)
  
}



#' Iterate through the grid regions of a spatial point pattern (starting from the one at the origin), finding 
#' the local PCFs within them. Ignores the marks of points. 
#' @param S spatial point pattern of class ppp
#' @param S_dim list of grid region breakpoints on S's grid: list(x_dim = breakpoints of the grid regions on x-axis, y_dim = breakpoints of the grid regions on y-axis). 
#' Assumes breakpoints are equally spaced for each of x_dim and y_dim, so the grid regions are equally-sized rectangles
#' @param bw bandwidth for pcf.ppp function
#' @param r vector of r-values for pcf.ppp function
#' @param min_cells if grid region has less than min_cells, then set PCF as missing
#' @return list of local window pcfs, one for each grid region
local_window_calc_pcfs <- function(S, S_dim, bw = NULL, r, min_cells = 2) {
  
  # iterate through the grid regions (starting from the one at the origin), finding 
  # the local pcf's within them
  
  ctr <- 1
  
  x_dim <- S_dim$x_dim
  y_dim <- S_dim$y_dim
  num_gp <- (length(x_dim) - 1) * (length(y_dim) - 1)
  
  loc_pcf_list <- vector(mode = "list", length = num_gp)
  
  local_area <- spatstat.geom::area(owin(xrange = c(x_dim[1], x_dim[2]),
                                         yrange = c(y_dim[1], y_dim[2])))
  
  for (j in 1:(length(x_dim) - 1)) { 
    
    for (i in 1:(length(y_dim) - 1)) {
      
      local_domain <- owin(xrange = c(x_dim[j], x_dim[j + 1]),
                           yrange = c(y_dim[i], y_dim[i + 1]))
      
      # clipping S onto the window local_domain, omitting all points outside it
      S_subset <- S
      Window(S_subset) <- local_domain
      
      # If local_domain is partially/fully outside window (which would make pcf.ppp display an error) or there are less than min_cells points, set pcf as missing
      if (!is.subset.owin(local_domain, Window(S)) | npoints(S_subset) < min_cells) { 
        
        loc_pcf_list[[ctr]] <- rep(NA, length(r))
        
      } else {
        
        flag <- F
        
        # if below returns a warning, it means that there are no interpoint distances <= max(r).
        # in this case, set the curve to missing. 
        tryCatch({res <- spatstat.explore::pcf.ppp(S_subset, correction = "translate",
                                                   ratio = FALSE, 
                                                   bw = bw, 
                                                   r = r)$trans}, warning = function(w) {flag <<- TRUE})
        if (flag) {
          loc_pcf_list[[ctr]] <- rep(NA, length(r))
        } else {
          loc_pcf_list[[ctr]] <- res
        }
        
      }
      
      ctr <- ctr + 1
      
    }
    
  }
  
  return(loc_pcf_list)
  
}



#' Square rooting, smoothing, and trimming the pcfs
#' @param r_grid is the set of r-values at which the pcfs are estimated.
#' @param trim_rvals To remove the pole at zero, functions will be trimmed so that only r-values greater than or equal to trim_rvals[1] will be included.
#' To localize the pcfs, functions will be trimmed so that only r-values less than or equal to trim_rvals[2] will be included
#' @param fv_obj_list the pcfs for all subjects
#' @return the processed pcfs for all subjects, as well as the corresponding trimmed r-sequence
smooth_local_pcfs <- function(r_grid, trim_rvals, fv_obj_list) {
  
  # removing the zero value, since it corresponds to Inf, which is not acceptable to smooth.spline
  r_grid <- r_grid[-1]
  
  N <- length(fv_obj_list)
  
  pcfs_smoothed <- vector(mode = "list", length = N)
  
  for (nn in 1:N) {
    
    pcfs <- fv_obj_list[[nn]][-1, ] # removing the Inf value, which is not acceptable to smooth.spline
    
    # square-rooting the functions to mitigate the skewness
    pcfs <- sqrt(pcfs)
    
    if (is.null(pcfs)) {
      
      pcfs_smoothed[[nn]] <- NULL
      
    } else {
      
      np <- ncol(pcfs)
      
      # estimate the function with smoothing.spline
      smooth_pcfs <- pcfs
      
      for (i in 1:np) {
        
        gc <- pcfs[, i]
        
        smooth_pcfs[, i] <- NA
        
        if (all(!is.na(gc))) {
          
          smooth_model <- smooth.spline(r_grid, gc, keep.data = F, cv = F)
          smooth_est <- smooth_model$y
          
          smooth_pcfs[, i] <- smooth_est
          
        }
        
      }
      
      # trimming the pcfs according to trim_rvals
      smooth_pcfs <- smooth_pcfs[r_grid >= trim_rvals[1] & r_grid <= trim_rvals[2], ]
      
      pcfs_smoothed[[nn]] <- smooth_pcfs
      
    }
    
  }
  
  
  
  # Trimming the r_grid according to trim_rvals
  r <- r_grid[r_grid >= trim_rvals[1] & r_grid <= trim_rvals[2]]
  
  return(list(pcfs_smoothed = pcfs_smoothed, r = r))
  
}
