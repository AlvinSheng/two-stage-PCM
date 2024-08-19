#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector clusterlabelupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                                 NumericVector Wtripletsum, const int nsites, NumericMatrix y_mat, NumericMatrix mu,
                                 int J, int M, double psi, NumericMatrix alpha_vec, IntegerVector cluster, NumericVector nu2, int g_idx)
{
  
  // Update the cluster labels
  // Create new objects
  int rowstart=0, rowend=0;
  double numnbr;
  double log_priorpottskernel, log_priory_matkernel;
  IntegerVector one2M;
  NumericVector log_probmasses(M), log_probmasses_diff(M), probmasses_ratio(M);
  double max_log_pm, log_py_matk_sum;
  NumericVector clusternew(nsites);
  
  one2M = seq_len(M);
  
  
  
  //  Update each cluster label in turn
  clusternew = cluster;
  for(int j = 0; j < nsites; j++)
  {
    
    // Find the probability masses for each possible cluster label m = 1, ..., M
    for(int m = 1; m <= M; m++) {
      
      // Evaluating kernel of Spatial Potts Model prior
      // needs count of neighbors with cluster label m
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      numnbr = 0;
      for(int l = rowstart; l < rowend; l++) numnbr += (clusternew[(Wtriplet(l,1) - 1)] == m);
      log_priorpottskernel = psi * numnbr + alpha_vec(m - 1, g_idx - 1);
      
      log_py_matk_sum = 0;
      for(int sj = 1; sj <= J; sj++) { 
        
        // Evaluating kernel of y_mat likelihood
        log_priory_matkernel = -1/(2*nu2[sj - 1]) * pow((y_mat(sj - 1, j) - mu(sj - 1, m - 1)), 2);
        
        log_py_matk_sum += log_priory_matkernel;
        
      }
      
      // Multiplying the kernels together
      log_probmasses[m - 1] = log_priorpottskernel + log_py_matk_sum;
      
    }
    
    // Use the log_probmasses to generate new cluster label
    
    max_log_pm = *std::max_element(log_probmasses.begin(), log_probmasses.end());
    
    log_probmasses_diff = log_probmasses;
    for(auto& element : log_probmasses_diff)
      element -= max_log_pm;
    
    probmasses_ratio = log_probmasses_diff;
    for(auto& element : probmasses_ratio)
      element = exp(element);
    
    // propose a value
    clusternew[j] = sample(one2M, 1, false, probmasses_ratio)[0];
    
  }
  
  return clusternew;
}



// [[Rcpp::export]]
NumericVector clusterlabelupdatesinglescore(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                                            NumericVector Wtripletsum, const int nsites, NumericVector y_vec, NumericVector mu,
                                            int J, int M, double psi, NumericMatrix alpha_vec, IntegerVector cluster, double nu2, int g_idx)
{
  
  // Update the cluster labels
  // Create new objects
  int rowstart=0, rowend=0;
  double numnbr;
  double log_priorpottskernel, log_priory_veckernel;
  IntegerVector one2M;
  NumericVector log_probmasses(M), log_probmasses_diff(M), probmasses_ratio(M);
  double max_log_pm;
  NumericVector clusternew(nsites);
  
  one2M = seq_len(M);
  
  
  
  //  Update each cluster label in turn
  clusternew = cluster;
  for(int j = 0; j < nsites; j++)
  {
    
    // Find the probability masses for each possible cluster label m = 1, ..., M
    for(int m = 1; m <= M; m++) {
      
      // Evaluating kernel of Spatial Potts Model prior
      // needs count of neighbors with cluster label m
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      numnbr = 0;
      for(int l = rowstart; l < rowend; l++) numnbr += (clusternew[(Wtriplet(l,1) - 1)] == m);
      log_priorpottskernel = psi * numnbr + alpha_vec(m - 1, g_idx - 1);
      
      // Evaluating kernel of y_vec likelihood
      log_priory_veckernel = -1/(2*nu2) * pow((y_vec[j] - mu[m - 1]), 2);
      
      // Multiplying the kernels together
      log_probmasses[m - 1] = log_priorpottskernel + log_priory_veckernel;
      
    }
    
    // Use the log_probmasses to generate new cluster label
    
    max_log_pm = *std::max_element(log_probmasses.begin(), log_probmasses.end());
    
    log_probmasses_diff = log_probmasses;
    for(auto& element : log_probmasses_diff)
      element -= max_log_pm;
    
    probmasses_ratio = log_probmasses_diff;
    for(auto& element : probmasses_ratio)
      element = exp(element);
    
    // propose a value
    clusternew[j] = sample(one2M, 1, false, probmasses_ratio)[0];
    
  }
  
  return clusternew;
}


