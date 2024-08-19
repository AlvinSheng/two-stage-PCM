# Monte Carlo Simulation Study

##### Generating Cluster Labels

# General Settings

M <- 3 # TRUE NUMBER OF CLUSTERS 

# TRUE BREAKPOINTS OF GRID REGIONS TO WHICH CLUSTERS WILL BE ASSIGNED
x_dim_true <- 0:12
y_dim_true <- 0:10

num_gp <- (length(x_dim_true) - 1) * (length(y_dim_true) - 1)



# Parameters of the Potts model

theta1 <- c(0, 0, 0, psi)
x <- matrix(1, nrow = length(y_dim_true) - 1, ncol = length(x_dim_true) - 1)
foo1 <- packPotts(x, M)





cluster_labs <- matrix(NA, nrow = num_subj, ncol = num_gp)

for (i in 1:num_subj) {
  
  out <- potts(foo1, theta1, nbatch = 10, blen = 20, boundary = "free")
  
  cluster_labs[i, ] <- c(unpackPotts(out$final))
  
}





##### Simulate Selected Spatial Point Patterns

# Calculate lambda2 from lambda1 and prob_mat, for the clusters. 
# First entry corresponds to the first intensity ratio, second entry corresponds to the second intensity ratio.
lambda2_vec <- c(lambda1 * (prob_mat[2, 1] / prob_mat[1, 1]), lambda1 * (prob_mat[2, 2] / prob_mat[1, 2]))

# Calculate unmarked intensities for the clusters
# First entry corresponds to the first intensity ratio, second entry corresponds to the second intensity ratio.
lambda_dot_vec <- c(lambda1 + lambda2_vec[1], lambda1 + lambda2_vec[2])

# Thomas parameters
thomas_kappa <- 5
thomas_scale <- 0.1

# Matérn I parameters (depends on lambda)
intens_root <- function(kappa, delta, intens) {kappa * exp(-kappa * pi * delta^2) - intens}
delta <- 0.037

maternI_kappa <- uniroot(intens_root, interval = c(0, 1), delta, intens = lambda_dot_vec[2], extendInt = "upX")$root



# Generating spatial point patterns 

ptm <- proc.time()

regime_spp <- vector(mode = "list", length = num_subj)

for (i in 1:num_subj) {
  
  cluster_im <- matrix(cluster_labs[i, ], nrow = length(y_dim_true) - 1, ncol = length(x_dim_true) - 1)
  
  if (1 %in% c(cluster_im)) {
    # CSR points
    w <- partial_window(mask = (cluster_im == 1))
    csr_ppp <- rpoispp(lambda = lambda_dot_vec[1], win = w)
    mark <- sample(c(1, 2), size = npoints(csr_ppp), replace = T, prob = prob_mat[, 1])
  } else {
    csr_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(csr_ppp) <- factor(mark, levels = c(1, 2))
  
  if (2 %in% c(cluster_im)) {
    # Thomas process points
    w <- partial_window(mask = (cluster_im == 2))
    thomas_ppp <- rThomas(kappa = thomas_kappa, scale = thomas_scale, mu = lambda_dot_vec[2] / thomas_kappa, win = w)
    mark <- sample(c(1, 2), size = npoints(thomas_ppp), replace = T, prob = prob_mat[, 2])
  } else{
    thomas_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(thomas_ppp) <- factor(mark, levels = c(1, 2))
  
  if (3 %in% c(cluster_im)) {
    # Matern I process points
    w <- partial_window(mask = (cluster_im == 3))
    maternI_ppp <- rMaternI(kappa = maternI_kappa, r = delta, win = w)
    mark <- sample(c(1, 2), size = npoints(maternI_ppp), replace = T, prob = prob_mat[, 2])
  } else{
    maternI_ppp <- ppp() # empty point pattern
    mark <- c()
  }
  marks(maternI_ppp) <- factor(mark, levels = c(1, 2))
  
  X <- superimpose(csr_ppp, thomas_ppp, maternI_ppp, 
                   W = owin(xrange = c(x_dim_true[1], tail(x_dim_true, 1)), yrange = c(y_dim_true[1], tail(y_dim_true, 1))))
  
  regime_spp[[i]] <- X
  
}

proc.time() - ptm




