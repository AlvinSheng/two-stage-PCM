library(here)
source(here("scripts/pkg_list.R"))



# Key Quantities 
M <- 3 # number of clusters
x_dim <- 12
y_dim <- 10
num_gp <- x_dim * y_dim

res_string <- paste0("res", y_dim, "x", x_dim)

print(paste0("Computations for ", res_string, ", M", M))



# Computation Strategy

# I will compute the mean estimate and the standard error estimate

# The foreach loop will work on these values of the psi_seq vector
psi_seq <- seq(0, 2.5, by = 0.01)
psi_range <- c(0, 2.5)
# and each combination of alpha values (alpha_2, ..., alpha_M). Each alpha falls within
alpha_range <- c(-5, 5)
# alpha_range is used for all alphas after the first; alpha_1 is set to zero in the potts model simulation
num_alpha_samp <- 10
# num_alpha_samp will be used to randomly select grid points to predict the canonical statistic for
num_realizations <- 100 # switch to 2 for quick runs
# number of Potts realizations used to estimate the canonical statistics
record_progress <- T
# whether to log progress in the intermediary_data/potts_canonical_stats folder
num_cores=detectCores() - 1 # strtoi(Sys.getenv(c("LSB_DJOB_NUMPROC")))
# number of cores to use



seed <- M * num_gp

iclust <- M - 1 

# sample nas grid points from multivariate normal to calculate canonical statistic for

nas <- num_alpha_samp*10^(iclust)

nas <- 1e5
# switch to 10 (or 1) for quick runs



cl <- makeCluster(num_cores) # not to overload your computer
registerDoParallel(cl)



if (record_progress) {
  writeLines(paste(paste0("g_alpha", 1:M, collapse = " "), "g_psi", paste0("c_alpha", 1:M, collapse = " "), "c_psi"), here("intermediary_data", "potts_canonical_stats", res_string, paste0("M", M, "ordered_grid_canon.txt")))
}

ptm <- proc.time()

# calculating spatial Potts samples
source(here("scripts/d_theta_scripts/d_theta_samps.R"), local = T)
stopCluster(cl)

proc.time() - ptm




