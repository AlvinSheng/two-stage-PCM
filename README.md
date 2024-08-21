# two-stage-PCM

Code for A Two-Stage Approach for Segmenting Spatial Point Patterns Applied to Tumor Immunology

by Alvin Sheng, Brian J Reich, Ana-Maria Staicu, Santhoshi N Krishnan, Arvind Rao, and Timothy L Frankel

We developed a method to segment the spatial point patterns (SPP) of cells in multiplex immunohistochemistry (mIF) images. The manuscript has been submitted to Biometrics.

The mIF images were obtained from patients at the University of Michigan Pancreatic Cancer Clinic who had undergone surgical resection for various pancreatic diseases.

## Key Scripts and Data

R Version 4.3.0 (Already Tomorrow) was used in writing and running the code.

### For evaluating the normalizing constant (running this step is optional; output is already generated)

Order of scripts to run:

1. scripts/d_theta_scripts/res10x12/canon_stat.R
3. scripts/d_theta_scripts/xgboost_potts_samps.R

Results in these intermediary data files:

1. intermediary_data/potts_canonical_stats/res10x12/M3ordered_grid_canon.txt
2. intermediary_data/potts_canonical_stats/res10x12/M3potts_samps_xgb.rds (already generated)

### Simulating spatial point patterns and applying the two-stage approach to it

Order of scripts to run:

1. scripts/pkg_list.R: ensure all of the packages in this script are installed and that this script can be run.
2. simulation_and_application_of_method.R
  * The script will load spatial_markov_model.R, which contains met_gibbs_potts() that runs the MCMC chains for the Potts clustering model (PCM).

Results in these intermediary data files:

1. intermediary_data/true_cluster_labs.rds: generated cluster labels which will be used to calculate the adjusted rand index (ARI)
2. intermediary_data/regime_spp.rds: generated SPPs
3. intermediary_data/pcfs_smoothed.rds: discretized vectors \hat{X}_{nl}

Results in this output:

1. MCMC_output/model_res.rds: an object with MCMC chains (or summaries of the chains) for all the parameters in the PCM, as well as model fitting information

### Running the PCM on the same simulated data, but assumming that the subjects are split into two groups

In simulation_and_application_of_method.R, change lines 90 and 101 according to the comments. 

Then, the function that runs the MCMC chain, met_gibbs_potts() in spatial_markov_model_group_idx.R, will assume by default that group = c(rep(1, floor(num_subj/2)), rep(2, ceiling(num_subj/2))), i.e., the first half of subjects are in group 1 and the second half are in group 2. Potts offsets will be computed separately for each group.


