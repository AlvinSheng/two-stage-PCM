# 11/22

sourceCpp(here("scripts/spatial_markov_model_group_idx.cpp"))



# Taken from duncanplee's code

common.Wcheckformat <- function(W)
{
  #### Check W is a matrix of the correct dimension
  if(!is.matrix(W)) stop("W is not a matrix.", call.=FALSE)
  n <- nrow(W)
  if(ncol(W)!= n) stop("W is not a square matrix.", call.=FALSE)    
  
  
  #### Check validity of inputed W matrix
  if(sum(is.na(W))>0) stop("W has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(W)) stop("W has non-numeric values.", call.=FALSE)
  if(min(W)<0) stop("W has negative elements.", call.=FALSE)
  if(sum(W!=t(W))>0) stop("W is not symmetric.", call.=FALSE)
  if(min(apply(W, 1, sum))==0) stop("W has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
  
  
  #### Create the triplet form
  ids <- which(W > 0, arr.ind = T)
  W.triplet <- cbind(ids, W[ids])
  W.triplet <- W.triplet[ ,c(2,1,3)]
  
  #W.triplet <- c(NA, NA, NA)
  #for(i in 1:n)
  #{
  #    for(j in 1:n)
  #    {
  #        if(W[i,j]>0)
  #        {
  #            W.triplet <- rbind(W.triplet, c(i,j, W[i,j]))     
  #        }else{}
  #    }
  #}
  #W.triplet <- W.triplet[-1, ]     
  n.triplet <- nrow(W.triplet) 
  W.triplet.sum <- tapply(W.triplet[ ,3], W.triplet[ ,1], sum)
  n.neighbours <- tapply(W.triplet[ ,3], W.triplet[ ,1], length)
  
  
  #### Create the start and finish points for W updating
  W.begfin <- cbind(c(1, cumsum(n.neighbours[-n])+1), cumsum(n.neighbours))
  #W.begfin <- array(NA, c(n, 2))     
  #temp <- 1
  #for(i in 1:n)
  #{
  #    W.begfin[i, ] <- c(temp, (temp + n.neighbours[i]-1))
  #    temp <- temp + n.neighbours[i]
  #}
  
  
  #### Return the critical quantities
  results <- list(W=W, W.triplet=W.triplet, n.triplet=n.triplet, W.triplet.sum=W.triplet.sum, n.neighbours=n.neighbours, W.begfin=W.begfin, n=n)
  return(results)   
}



# Taken from duncanplee's code
# Acceptance rates - maximum limit on the proposal sd
common.accceptrates2 <- function(accept, sd, min, max, sd.max)
{
  #### Update the proposal standard deviations
  rate <- 100 * accept[1] / accept[2]
  
  if(rate > max)
  {
    sd <- sd + 0.2 * sd
    sd[which(sd>sd.max)] <- sd.max
  }else if(rate < min)              
  {
    sd <- sd - 0.2 * sd
  }else
  {
  }
  
  return(sd)
}



# Taken from duncanplee's code
# Compute the DIC. WAIC,LMPL and loglikelihood
common.modelfit <- function(samples.loglike, deviance.fitted)
{
  #### WAIC
  p.w <- sum(apply(samples.loglike, 2, var), na.rm=TRUE)
  mean.like <- apply(exp(samples.loglike),2,mean)
  mean.min <- min(mean.like[mean.like>0])
  mean.like[mean.like==0] <- mean.min
  lppd <- sum(log(mean.like), na.rm=TRUE)
  WAIC <- -2 * (lppd - p.w)
  
  
  #### Compute the Conditional Predictive Ordinate
  CPO <- 1/apply(exp(-samples.loglike), 2, mean)
  mean.min <- min(CPO[CPO>0])
  CPO[CPO==0] <- mean.min
  LMPL <- sum(log(CPO), na.rm=TRUE)    
  
  
  #### DIC
  mean.deviance <- -2 * sum(samples.loglike, na.rm=TRUE) / nrow(samples.loglike)
  p.d <- mean.deviance - deviance.fitted
  DIC <- deviance.fitted + 2 * p.d
  
  
  #### loglikelihood
  loglike <- -0.5 * deviance.fitted
  
  
  #### Model fit criteria
  modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike)
  names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")
  return(modelfit)  
}



# Taken from duncanplee's code
# Compute the DIC. WAIC,LMPL and loglikelihood
common.modelfit.summarized.loglike <- function(mean.like, mean.recip.like, var.loglike, sum.loglike, n.keep, deviance.fitted)
{
  #### WAIC
  p.w <- sum(var.loglike, na.rm=TRUE)
  mean.min <- min(mean.like[mean.like>0])
  mean.like[mean.like==0] <- mean.min
  lppd <- sum(log(mean.like), na.rm=TRUE)
  WAIC <- -2 * (lppd - p.w)
  
  
  #### Compute the Conditional Predictive Ordinate
  CPO <- 1/mean.recip.like
  mean.min <- min(CPO[CPO>0])
  CPO[CPO==0] <- mean.min
  LMPL <- sum(log(CPO), na.rm=TRUE)    
  
  
  #### DIC
  mean.deviance <- -2 * sum.loglike / n.keep
  p.d <- mean.deviance - deviance.fitted
  DIC <- deviance.fitted + 2 * p.d
  
  
  #### loglikelihood
  loglike <- -0.5 * deviance.fitted
  
  
  #### Model fit criteria
  modelfit <- c(DIC, p.d, WAIC, p.w, LMPL, loglike)
  names(modelfit) <- c("DIC", "p.d", "WAIC", "p.w", "LMPL", "loglikelihood")
  return(modelfit)  
}



# similar to quadform, but phi is replaced by cluster labels
quad_mat <- function(M, cluster) {
  
  qmat <- matrix(0, nrow = M, ncol = M)
  
  for (m in 1:M) {
    qmat[m, m] <- sum(cluster == m)
  }
  
  return(qmat)
  
}



# similar to quad_mat, but involves both the cluster labels and y_vec
quad_vec <- function(M, cluster, y_vec) {
  
  qvec <- rep(NA, length = M) 
  
  for (m in 1:M) {
    qvec[m] <- sum(y_vec[cluster == m])
  }
  
  return(qvec)
  
} 



# Turning vector into dummy-coded matrix. Not strictly necessary, just for
# interfacing with duncanplee's code.
# Assumes clusters are numbered from 1, ..., M. 
cluster2dummy <- function(M, cluster) {
  
  K <- length(cluster)
  
  Xc <- matrix(0, nrow = K, ncol = M)
  
  for (k in 1:K) {
    
    Xc[k, cluster[k]] <- 1
    
  }
  
  return(Xc)
  
}





#' Metropolis-Gibbs algorithm for Potts clustering model
#'
#' @param y_array N by J by num_gp matrix of response values (J is the number of PC scores)
#' @param M number of possible clusters
#' @param W adjacency matrix (assuming no additional islands)
#' @param n_burn_in number of burn-in iterations
#' @param n_iter number of kept iterations
#' @param thin level of thinning to apply to the MCMC samples
#' @param keep_only vector of subjects to keep all parameters for; for other subjects, 
#' keep track of the running posterior mean. If length(keep_only) < number of subjects, 
#' then cluster_proportion_array is also returned.
#' @param potts_samps object needed to compute the log(d(theta)) value for the given M value, where d(theta) is the normalizing constant of the spatial potts model.
#' First item in list is the model for predicting the psi canonical statistic curve, and the second item in list are models for predicting the alpha canonical statistic curve.  
#' All functions assume that the alpha parameters are in order from least to greatest. 
#' Assumes that the functions correspond to the resolution of the adjacency matrix W
#' @param true_cluster_labs for the purpose of the simulation studies, the adjusted Rand index 
#' will be calculated within the MCMC loop to measure performance, while adjusting for label switching
#' @param finer_dim list of breakpoints (x_dim, y_dim) for true_cluster_labs if it was generated at a finer resolution than the model assumes
#' @param label_unswitch reorder the labels within the MCMC loop based on the Pivotal Reordering Algorithm proposed by Marin, Mengersen, and Robert
#' @param mu_inform (optional) this is the prior precision for mu. Set this to a big value to make a more informative Normal prior for mu centered at the cluster means from k-means. If NULL, use a noninformative Normal prior with the sample mean and sample variance.
#' @param psi_fix (optional) fix the psi parameter at this value
#' @param alpha_fix if NULL, let alphas be estimated via MCMC. If vector given, freeze alphas to that.
#' @param group factor indicating the group(s) the subjects are in, whose alpha parameters will be estimated separately
#' @return MCMC chain.
#' Notes on the return list.
#' mean.y.mat: this is the given y_array, but the missing values are replaced by the mean of the posterior predictive distribution
#' mean.resid: this is the posterior mean of the difference between y_array (missing values imputed) and fitted value mu 
#'
#' Example usage.
#' set.seed(101)
#' my_res_pcfs <- met_gibbs_potts(y_array=y_array, M=3, W=W,
#'                                n_burn_in = 10000, n_iter = 20000, thin = 1,
#'                                keep_only = c(1, 31), potts_samps = potts_samps, true_cluster_labs = cluster_labs, finer_dim = finer_dim,
#'                                label_unswitch = T, mu_inform = NULL, psi_fix = NULL, alpha_fix = NULL, group = factor(c(rep("grp1", floor(dim(y_array)[1]/2)), rep("grp2", ceiling(dim(y_array)[1]/2)))))
#' saveRDS(my_res_pcfs, file = here("my_res_pcfs.rds"))
#' @export
met_gibbs_potts <- function(y_array, M, W, n_burn_in, n_iter, thin = 1, 
                            keep_only = 1:dim(y_array)[1], potts_samps, true_cluster_labs = NULL, finer_dim = NULL, label_unswitch = T, mu_inform = NULL, psi_fix = NULL, alpha_fix = NULL, 
                            group = c(rep(1, floor(dim(y_array)[1]/2)), rep(2, ceiling(dim(y_array)[1]/2)))) {
  
  N <- dim(y_array)[1]
  J <- dim(y_array)[2]
  K <- dim(y_array)[3]
  S <- n_iter
  
  # Initial values
  ## for nu2
  nu2 <- apply(y_array, 2, function(mat) var(c(mat), na.rm = T))
  
  ## for y_array
  which.na <- is.na(y_array)
  n.na <- sum(which.na)
  ## mean imputation of missing Y values
  y_array_da <- y_array
  resid <- y_array # will be updated within for-loop, right before being stored in the MCMC array
  
  for (sn in 1:dim(y_array)[1]) {
    for (sj in 1:dim(y_array)[2]) {
      vec_isna <- is.na(y_array_da[sn, sj, ])
      y_array_da[sn, sj, ][vec_isna] <- mean(y_array_da[sn, sj, ], na.rm = T)
    }
  }
  
  ## for cluster 
  y_array_flatten <- apply(y_array_da, 2, c)
  kmeans_res <- kmeans(y_array_flatten, M)
  mu.mean <- kmeans_res$centers
  cluster_mat <- matrix(kmeans_res$cluster, nrow = N, ncol = K)
  
  ## for the Potts parameters
  ### for psi
  psi_seq <- attr(potts_samps, "psi_seq", exact = T)
  if(is.null(psi_seq)) {psi_seq <- c(seq(0, 2.4, by = 0.2), 2.5)}
  min_psi <- min(psi_seq)
  max_psi <- max(psi_seq)
  
  psi <- 0 # initial psi
  
  ### for psi spline calculations for log(d(theta(psi)))
  p <- attr(potts_samps, "p", exact = T)
  if(is.null(p)) {p <- 10}
  
  xgb_mod_psis <- potts_samps[[1]]
  
  PSI <- matrix(psi_seq, length(psi_seq), p, byrow = F) 
  pow <- matrix(1:p, length(psi_seq), p, byrow = T)
  X <- PSI^pow
  
  if (!is.null(psi_fix)) {
    psi <- psi_fix
  }
  
  ### for alphas
  alpha_seq <- attr(potts_samps, "alpha_seq", exact = T)
  if(is.null(alpha_seq)) {alpha_seq <- seq(-5, 5, by = 0.05)}
  min_alpha <- min(alpha_seq)
  max_alpha <- max(alpha_seq)
  
  xgb_mod_alphas <- potts_samps[[2]]
  
  group_idx <- as.numeric(group)
  num_group <- length(unique(group))
  
  if (!is.null(alpha_fix)) {
    alpha_vec <- alpha_fix
  } else {
    alpha_vec <- do.call("cbind", lapply(1:num_group, function(x) rep(0, M))) # rep(0, M) # initial alpha1, ..., alphaM
  }
  
  ## for mu 
  Xc <- cluster2dummy(M, cluster = c(cluster_mat))
  mu <- matrix(NA, nrow = J, ncol = M)
  for (sj in 1:J) {
    mu[sj, ] <- mu.mean[, sj]
  }
  
  #### Matrices to store samples
  n.keep <- floor(n_iter / thin)
  samples.mu <- array(NA, c(n.keep, J, M))
  samples.cluster <- array(NA, c(n.keep, length(keep_only), K))
  samples.nu2 <- array(NA, c(n.keep, J))
  samples.psi <- array(NA, c(n.keep))
  samples.alpha <- array(NA, c(n.keep, M-1, num_group)) # alpha2, ..., alphaM
  samples.loglike <- array(NA, c(n.keep, length(keep_only), J, K))
  
  cluster_proportion_array <- array(0, c(N, K, M))
  mean.y.mat <- array(0, c(N, J, K))
  mean.resid <- array(0, c(N, J, K)) # Effective sizes/traceplots for residuals are not needed, so tracking the mean residuals is good enough
  if (!is.null(true_cluster_labs)) {
    mean.ARI <- 0
  }
  
  # If length(keep_only) is less than N, then track the average of parameters across subjects
  # or other summary statistics, in the case of log-likelihood
  if (length(keep_only) < N) {
    mean.loglike <- array(0, c(N, J, K))
    second.mom.loglike <- array(0, c(N, J, K))
    mean.like <- array(0, c(N, J, K))
    mean.recip.like <- array(0, c(N, J, K))
  }
  
  #### adjacency matrix quantities
  W.quants <- common.Wcheckformat(W)
  W <- W.quants$W
  W.triplet <- W.quants$W.triplet
  n.triplet <- W.quants$n.triplet
  W.triplet.sum <- W.quants$W.triplet.sum
  n.neighbours <- W.quants$n.neighbours 
  W.begfin <- W.quants$W.begfin
  
  # find x_dim_coarse, y_dim_coarse from W if true_cluster_labs was 
  # generated at a finer resolution than the model assumes
  if (!is.null(finer_dim)) {
    # assumes W corresponds to the adjacencies of a rectangular grid in column-major order
    y_len_coarse <- which(W[, 1] == 1)[2] - 1
    x_len_coarse <- nrow(W) / y_len_coarse
    y_dim_coarse <- seq(finer_dim$y_dim[1], tail(finer_dim$y_dim, 1), length.out = y_len_coarse + 1)
    x_dim_coarse <- seq(finer_dim$x_dim[1], tail(finer_dim$x_dim, 1), length.out = x_len_coarse + 1)
  }
  
  # Prior quantities
  ## mu update quantities (Gibbs)
  if (is.null(mu_inform)) {
    bmu <- 1 / apply(y_array, 2, function(mat) var(c(mat), na.rm = T))
    precision_diags <- lapply(bmu, function(x) diag(rep(x, M)))
    prior.precision.mu <- abind(precision_diags, along = 3)
  } else {
    prior.precision.mu <- array(diag(rep(mu_inform, M)), dim = c(M, M, J))
  }
  if (is.null(mu_inform)) {
    amu <- apply(y_array, 2, function(mat) mean(c(mat), na.rm = T))
    mean_reps <- lapply(amu, function(x) rep(x, M)) # matrix(mean(c(y_array), na.rm = T), J, M) 
    prior.mean.mu <- do.call("rbind", mean_reps)
  } else {
    prior.mean.mu <- t(mu.mean)
  }
  
  ## nu2 update quantities (Gibbs)
  prior.nu2 <- c(1, 0.01) # prior shape and scale
  nu2.posterior.shape <- prior.nu2[1] + N * K / 2
  ## potts update quantities (Metropolis-Hastings)
  accept_psi <- rep(0,2)
  proposal.sd.psi <- 2
  accept_alpha <- array(0, c(2, M-1, num_group))
  proposal.sd.alphas <- matrix(4, nrow = M-1, ncol = num_group)
  ## for simulated annealing
  temperature <- 1 / n_burn_in 
  
  ## other preliminary values
  num_nbr_mat <- diag(W.quants$n.neighbours)
  
  mu_array <- array(NA, dim = dim(y_array)) 
  for (sj in 1:J) {
    mu_array[, sj, ] <- mu.mean[, sj][cluster_mat]
  }
  ## for PRA
  kmeans_likel <- 0
  for (sj in 1:J) {
    kmeans_likel <- kmeans_likel + sum(dnorm(y_array[, sj, ], mean = mu_array[, sj, ], sd = sqrt(nu2[sj]), log = TRUE), na.rm=TRUE)
  }
  if (label_unswitch) {
    init_pivot <- t(mu.mean)
    pivot <- init_pivot
    best_likel <- kmeans_likel
  }
  
  # calculate cluster_count for the M groups to be freezed during burn-in
  pairwise_count_init <- 0
  for (sn in 1:N) {
    pairwise_count_init <- pairwise_count_init + pairwise_same_cluster_count(cluster_mat[sn, ], W.triplet)
  }
  cluster_count_init <- matrix(NA, nrow = M, ncol = num_group)
  for (sg in 1:num_group) {
    for (sm in 1:M) {
      ### current cluster count, canonical statistic
      cluster_count_init[sm, sg] <- sum(cluster_mat[group_idx == sg, ] == sm)
    }
  }
  
  
  
  for (s in 1:(S + n_burn_in)) {
    # Gibbs/Metropolis for mu, nu2, cluster, psi
    ## mu
    quad_mat_sum <- 0
    for (sn in 1:N) {
      quad_mat_sum <- quad_mat_sum + quad_mat(M, cluster_mat[sn, ])
    }
    for (sj in 1:J) {
      fc.precision <- prior.precision.mu[, , sj] + quad_mat_sum / nu2[sj]
      fc.var <- solve(fc.precision)
      quad_vec_sum <- 0
      for (sn in 1:N) {
        quad_vec_sum <- quad_vec_sum + quad_vec(M, cluster_mat[sn, ], y_array_da[sn, sj, ])
      }
      mu.offset2 <- quad_vec_sum / nu2[sj] + prior.precision.mu[, , sj] %*% prior.mean.mu[sj, ]
      fc.mean <- fc.var %*% mu.offset2
      chol.var <- t(chol(fc.var))
      mu[sj, ] <- fc.mean + chol.var %*% rnorm(M)
    }
    
    ## cluster 
    for (sn in 1:N) {
      if (J > 1) {
        cluster_mat[sn, ] <- clusterlabelupdate(W.triplet, W.begfin, W.triplet.sum, nsites = K, y_mat = y_array_da[sn, , ],
                                                mu = mu, J, M, psi, alpha_vec, cluster_mat[sn, ], nu2, g_idx = group_idx[sn])
      } else {
        cluster_mat[sn, ] <- clusterlabelupdatesinglescore(W.triplet, W.begfin, W.triplet.sum, nsites = K, y_array_da[sn, , ],
                                                           mu = mu, J, M, psi, alpha_vec, cluster_mat[sn, ], nu2, g_idx = group_idx[sn])
      }
    }
    # to make the cluster memberships identifiable, relabel the clusters based on PRA:
    # find the permutation that maximizes the dot product between mu and the pivot
    if (label_unswitch) {
      if (s > n_burn_in) {
        current_pivot <- pivot
      } else {
        current_pivot <- init_pivot
      }
      perms_mat <- permutations(n = M, r = M) # iterate through all the possible permutations
      best_perm <- 1:M
      best_val <- -Inf
      for (i in 1:nrow(perms_mat)) {
        perm <- perms_mat[i, ]
        dotp <- sum(mu[, perm] * current_pivot)
        if(dotp >= best_val) {
          best_val <- dotp
          best_perm <- perm
        }
      }
      # using the permutation that maximized the dot product
      cluster_mat <- apply(cluster_mat, c(1, 2), function(ele) best_perm[ele])
    }
    
    if (!is.null(true_cluster_labs)) {
      if (!is.null(finer_dim)) { # if true_cluster_labs was generated at a finer resolution
        mcmc_ARI <- adjustedRandIndex(cluster_labs_finer(cluster_mat, x_dim_coarse, y_dim_coarse, finer_dim), 
                                      true_cluster_labs)
      } else {
        mcmc_ARI <- adjustedRandIndex(cluster_mat, true_cluster_labs)
      }
    }
    
    ## psi
    if (is.null(psi_fix)) {
      ### current pairwise same cluster count, canonical statistic
      if (s <= n_burn_in) {
        pairwise_count <- pairwise_count_init
      } else {
        pairwise_count <- 0
        for (sn in 1:N) {
          pairwise_count <- pairwise_count + pairwise_same_cluster_count(cluster_mat[sn, ], W.triplet)
        }
      }
      
      # assembling matrix for prediction
      # in the order of other alphas (if any) and psi
      psi_seq_mat <- cbind(matrix(alpha_vec[, 1], nrow = length(psi_seq), ncol = M, byrow = T), psi_seq)
      preds_psi <- rep(NA, length(psi_seq))
      offsets_order <- c(order(alpha_vec[, 1]), M+1)
      # using the models based on the rank of the alphas
      for (i in 1:length(psi_seq)) {
        pt_ordered <- psi_seq_mat[i, offsets_order, drop = F]
        colnames(pt_ordered) <- c(paste0("alpha", 1:M), "psi")
        
        xgb_test = xgb.DMatrix(data = pt_ordered)
        preds_psi[i] <- predict(xgb_mod_psis, xgb_test)
      }
      beta1 <- lm(preds_psi ~ X)$coef
      
      # doing it again for group 2
      psi_seq_mat <- cbind(matrix(alpha_vec[, 2], nrow = length(psi_seq), ncol = M, byrow = T), psi_seq)
      preds_psi <- rep(NA, length(psi_seq))
      offsets_order <- c(order(alpha_vec[, 2]), M+1)
      # using the models based on the rank of the alphas
      for (i in 1:length(psi_seq)) {
        pt_ordered <- psi_seq_mat[i, offsets_order, drop = F]
        colnames(pt_ordered) <- c(paste0("alpha", 1:M), "psi")
        
        xgb_test = xgb.DMatrix(data = pt_ordered)
        preds_psi[i] <- predict(xgb_mod_psis, xgb_test)
      }
      beta2 <- lm(preds_psi ~ X)$coef
      
      # interpolation
      X1   <- c(psi, (psi^(1:p+1))/(1:p+1))
      logd_current <- - sum(group_idx == 1) * sum(X1*beta1) - sum(group_idx == 2) * sum(X1*beta2)
      logprob_psi_current <- logd_current + psi * pairwise_count
      
      psi_proposal <- rtruncnorm(n=1, a=min_psi, b=max_psi, mean=psi, sd=proposal.sd.psi)
      X2   <- c(psi_proposal, (psi_proposal^(1:p+1))/(1:p+1))
      logd_proposal <- - sum(group_idx == 1) * sum(X2*beta1) - sum(group_idx == 2) * sum(X2*beta2)
      logprob_psi_proposal <- logd_proposal + psi_proposal * pairwise_count
      hastings <- log(dtruncnorm(x=psi, a=min_psi, b=max_psi, mean=psi_proposal, sd=proposal.sd.psi)) - log(dtruncnorm(x=psi_proposal, a=min_psi, b=max_psi, mean=psi, sd=proposal.sd.psi))
      prob <- exp(logprob_psi_proposal - logprob_psi_current + hastings)
      # simulated annealing
      prob <- prob^temperature
      ## Accept or reject the proposal
      if(prob > runif(1))
      {
        psi <- psi_proposal
        accept_psi[1] <- accept_psi[1] + 1
      }
      accept_psi[2] <- accept_psi[2] + 1
    }
    
    ## alphas
    if (!is.null(alpha_fix)) {
      alpha_vec <- alpha_fix
    } else {
      for (sm in 2:M) {
        for (sg in 1:num_group) {
          ### current cluster count, canonical statistic
          cluster_count <- sum(cluster_mat[group_idx == sg, ] == sm)
          if (s <= n_burn_in) {
            cluster_count <- cluster_count_init[sm, sg] # Stabilizing the alpha estimates by freezing cluster counts at those calculated by k-means.
          }
          logprob_alpha_current <- alpha_vec[sm, sg] * cluster_count - alpha_vec[sm, sg]^2/(2*.25)
          alpha_proposal <- rtruncnorm(n=1, a=min_alpha, b=max_alpha, mean=alpha_vec[sm, sg], sd=proposal.sd.alphas[sm-1, sg])
          logprob_alpha_proposal <- alpha_proposal * cluster_count - alpha_proposal^2/(2*.25)
          
          # assembling matrix for prediction
          # combining alpha_seq with other_parms
          alpha_bw <- seq(alpha_vec[sm, sg], alpha_proposal, length = length(psi_seq))
          alpha_bw_mat <- cbind(matrix(alpha_vec[, sg], nrow = length(alpha_bw), ncol = M, byrow = T), rep(psi, length(alpha_bw)))
          alpha_bw_mat[, sm] <- alpha_bw
          preds_alpha <- rep(NA, length(alpha_bw))
          # using the models based on the rank of the alphas
          for (i in 1:length(alpha_bw)) {
            pt <- as.numeric(alpha_bw_mat[i, ])
            offsets_order <- c(order(pt[1:M]), M+1)
            pt_ordered <- alpha_bw_mat[i, offsets_order, drop = F]
            colnames(pt_ordered) <- c(paste0("alpha", 1:M), "psi")
            
            xgb_test = xgb.DMatrix(data = pt_ordered)
            preds_alpha[i] <- predict(xgb_mod_alphas[[which(offsets_order == sm)]], xgb_test)
            if (i > 1) {
              if (preds_alpha[i] < preds_alpha[i - 1]) {
                preds_alpha[i] <- preds_alpha[i - 1]
              }
            }
          }
          logdiff <- piecewise(alpha_bw, preds_alpha)$area
          
          hastings <- log(dtruncnorm(x=alpha_vec[sm, sg], a=min_alpha, b=max_alpha, mean=alpha_proposal, sd=proposal.sd.alphas[sm-1, sg])) - log(dtruncnorm(x=alpha_proposal, a=min_alpha, b=max_alpha, mean=alpha_vec[sm, sg], sd=proposal.sd.alphas[sm-1, sg]))
          prob <- exp(- sum(group_idx == sg) * logdiff + logprob_alpha_proposal - logprob_alpha_current + hastings)
          # simulated annealing
          prob <- prob^temperature
          ## Accept or reject the proposal
          if(prob > runif(1))
          {
            alpha_vec[sm, sg] <- alpha_proposal
            accept_alpha[1, sm-1, sg] <- accept_alpha[1, sm-1, sg] + 1
          }
          accept_alpha[2, sm-1, sg] <- accept_alpha[2, sm-1, sg] + 1
        }
      }
    }
    
    ## nu2
    for (sj in 1:J) { # re-used in data augmentation step below
      mu_array[, sj, ] <- mu[sj, ][cluster_mat]
    }
    nu2.posterior.scale <- prior.nu2[2] + 0.5 * apply(y_array_da - mu_array, 2, function(mat) {sum(mat^2)})
    nu2 <- sapply(nu2.posterior.scale, function(x) {1 / rgamma(1, nu2.posterior.shape, scale = (1 / x))})
    
    # Sampling for missing y_list values - data augmentation
    if (n.na > 0) {
      for(sj in 1:J){
        mat_isna <- is.na(y_array[, sj, ]) 
        y_array_da[, sj, ][mat_isna] <- rnorm(n = sum(mat_isna),
                                              mean = mu_array[, sj, ][mat_isna],
                                              sd = sqrt(nu2[sj]))
        resid[, sj, ] <- y_array_da[, sj, ] - mu_array[, sj, ]
      }
    }
    ## calculate the deviance 
    loglike <- array(NA, c(N, J, K))
    for (sj in 1:J) {
      loglike[, sj, ] <- dnorm(y_array[, sj, ], mean = mu_array[, sj, ], sd = sqrt(nu2[sj]), log = TRUE)
    }
    
    if (s > n_burn_in)
    {
      if ((s - n_burn_in) %% thin == 0) {
        ele <- (s - n_burn_in) / thin
        samples.mu[ele, , ] <- mu
        samples.cluster[ele, , ] <- cluster_mat[keep_only, ]
        samples.nu2[ele, ] <- nu2
        samples.psi[ele] <- psi
        samples.alpha[ele, , ] <- alpha_vec[-1, ]
        samples.loglike[ele, , , ] <- loglike[keep_only, , ]
      }
      for (sm in 1:M) {
        cluster_proportion_array[, , sm] <- cluster_proportion_array[, , sm] + (cluster_mat == sm) / n_iter
      }
      mean.y.mat <- mean.y.mat + y_array_da / n_iter
      mean.resid <- mean.resid + resid / n_iter
      if (!is.null(true_cluster_labs)) {
        mean.ARI <- mean.ARI + mcmc_ARI / n_iter
      }
      if (length(keep_only) < N) {
        mean.loglike <- mean.loglike + loglike / n_iter
        second.mom.loglike <- second.mom.loglike + loglike^2 / n_iter
        mean.like <- mean.like + exp(loglike) / n_iter
        mean.recip.like <- mean.recip.like + exp(-loglike) / n_iter
      }
    }
    
    # simulated annealing
    if (s <= n_burn_in) {
      temperature <- s/n_burn_in
    }
    
    # Self-tune the acceptance probabilities
    if(ceiling(s/100) == floor(s/100) & s <= n_burn_in)
    {
      if (is.null(psi_fix)) {
        proposal.sd.psi <- common.accceptrates2(accept_psi[1:2], proposal.sd.psi, 30, 50, sd.max = 2.5)
        accept_psi <- c(0,0)
      }
      if (is.null(alpha_fix)) {
        for (sg in 1:num_group) {
          for (sm in 2:M) {
            proposal.sd.alphas[sm-1, sg] <- common.accceptrates2(accept_alpha[, sm-1, sg], proposal.sd.alphas[sm-1, sg], 30, 50, sd.max = 5)
            accept_alpha[, sm-1, sg] <- c(0,0)
          }
        }
      }
    }
    
    if (s <= n_burn_in) {
      # Find the pivot corresponding to the maximized likelihood during the burn-in phase
      burn_in_likel <- sum(loglike, na.rm = T)
      if (label_unswitch) {
        if (burn_in_likel >= best_likel) {
          pivot <- mu
          best_likel <- burn_in_likel
        }
      }
    }
    
  }
  
  #### Compute the acceptance rates
  accept_psi_final <- 100 * accept_psi[1] / accept_psi[2]
  accept_alpha_final <- c(100 * accept_alpha[1, , 1] / accept_alpha[2, , 1], 100 * accept_alpha[1, , 2] / accept_alpha[2, , 2])
  accept_final <- c(accept_psi_final, accept_alpha_final)
  names(accept_final) <- c("psi", paste0("alpha_firstgrp", 2:M), paste0("alpha_secondgrp", 2:M))
  
  #### Compute the fitted deviance
  nu2.mean <- apply(samples.nu2, 2, function(vec) mean(vec, na.rm = T))
  mu.mean.post <- apply(samples.mu, c(2, 3), function(vec) mean(vec, na.rm = T))
  mode_cluster_mat <- apply(cluster_proportion_array, c(1, 2), which.max)
  mu.mean.array <- array(NA, dim = dim(y_array_da)) # re-used in data augmentation step below
  for (sj in 1:J) {
    mu.mean.array[, sj, ] <- mu.mean.post[sj, ][mode_cluster_mat]
  }
  deviance_vec <- rep(NA, J)
  for (sj in 1:J) {
    deviance_vec[sj] <- -2 * sum(dnorm(y_array[, sj, ], mean = mu.mean.array[, sj, ], sd = sqrt(nu2.mean[sj]), log = TRUE), na.rm=TRUE)
  }
  deviance.fitted <- sum(deviance_vec)
  
  #### Model fit criteria
  # flatten samples.loglike
  if (length(keep_only) == N) { 
    samples.loglike.flatten.list <- vector(mode = "list", length = J)
    for (sj in 1:J) {
      samples.loglike.flatten <- matrix(NA, nrow = n.keep, ncol = K * N) 
      for (sn in 1:N) {
        samples.loglike.flatten[, ((sn - 1) * K + 1):(sn * K)] <- samples.loglike[, sn, sj, ]
      }
      samples.loglike.flatten.list[[sj]] <- samples.loglike.flatten
    }
    samples.loglike.flatten.NJ <- do.call("cbind", samples.loglike.flatten.list)
    modelfit <- common.modelfit(samples.loglike.flatten.NJ, deviance.fitted) 
  } else {
    var.loglike <- (n.keep / (n.keep - 1)) * (second.mom.loglike - mean.loglike^2)
    sum.loglike <- sum(mean.loglike * n.keep, na.rm = T)
    modelfit <- common.modelfit.summarized.loglike(mean.like = mean.like, mean.recip.like = mean.recip.like, 
                                                   var.loglike = var.loglike, sum.loglike = sum.loglike, 
                                                   n.keep = n.keep, deviance.fitted = deviance.fitted) 
  }
  
  if (length(keep_only) == N) {
    samples <- list(mu = samples.mu, cluster_mat = samples.cluster, 
                    nu2 = mcmc(samples.nu2), psi = mcmc(samples.psi), alphas = samples.alpha)
  } else {
    samples <- list(mu = samples.mu, cluster_mat = samples.cluster, 
                    cluster_proportion_array = cluster_proportion_array, nu2 = mcmc(samples.nu2), psi = mcmc(samples.psi), alphas = samples.alpha)
  }
  
  if (!is.null(true_cluster_labs)) {
    results <- list(samples = samples, mean.y.mat = mean.y.mat, mean.resid = mean.resid, modelfit = modelfit, mean.ARI = mean.ARI, accept = accept_final)
  } else {
    results <- list(samples = samples, mean.y.mat = mean.y.mat, mean.resid = mean.resid, modelfit = modelfit, accept = accept_final)
  }
  
  return(results)
}



