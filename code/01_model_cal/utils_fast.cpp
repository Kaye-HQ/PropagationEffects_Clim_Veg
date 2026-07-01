#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector roll_mean_right_partial_cpp(NumericVector x, int window) {
  int n = x.size();
  NumericVector res(n, NA_REAL);
  
  for(int i = 0; i < n; i++) {
    double sum = 0.0;
    int count = 0;
    // Right-aligned window
    int start_idx = std::max(0, i - window + 1);
    
    for(int j = start_idx; j <= i; j++) {
      if(!NumericVector::is_na(x[j])) {
        sum += x[j];
        count++;
      }
    }
    if(count > 0) res[i] = sum / count;
  }
  return res;
}

// [[Rcpp::export]]
List methods_evaluate_cpp(String method, NumericVector obs, NumericMatrix sim) {
  int n_sims = sim.ncol();
  int n_obs = obs.size();
  NumericVector value(n_sims, NA_REAL), fit(n_sims, NA_REAL);
  
  double mean_obs = 0, sd_obs = 0;
  int count_obs = 0;
  for(int i=0; i<n_obs; i++) {
    if(!NumericVector::is_na(obs[i])) { mean_obs += obs[i]; count_obs++; }
  }
  mean_obs /= count_obs;
  
  for(int i=0; i<n_obs; i++) {
    if(!NumericVector::is_na(obs[i])) { sd_obs += pow(obs[i] - mean_obs, 2); }
  }
  sd_obs = sqrt(sd_obs / (count_obs - 1));
  double SS_total = sd_obs * sd_obs * (count_obs - 1);
  
  for(int j=0; j<n_sims; j++) {
    NumericVector cur_sim = sim(_, j);
    double mean_sim = 0, sd_sim = 0;
    int count_sim = 0;
    
    for(int i=0; i<n_obs; i++) {
      if(!NumericVector::is_na(cur_sim[i])) { mean_sim += cur_sim[i]; count_sim++; }
    }
    mean_sim /= count_sim;
    
    for(int i=0; i<n_obs; i++) {
      if(!NumericVector::is_na(cur_sim[i])) { sd_sim += pow(cur_sim[i] - mean_sim, 2); }
    }
    sd_sim = sqrt(sd_sim / (count_sim - 1));
    
    if (method == "KGE") {
      double cov = 0;
      for(int i=0; i<n_obs; i++) {
        if(!NumericVector::is_na(obs[i]) && !NumericVector::is_na(cur_sim[i])) {
          cov += (obs[i] - mean_obs) * (cur_sim[i] - mean_sim);
        }
      }
      double r = cov / ((count_obs - 1) * sd_obs * sd_sim);
      double alpha = sd_sim / sd_obs;
      double beta = mean_sim / mean_obs;
      
      value[j] = 1.0 - sqrt(pow(r - 1.0, 2) + pow(alpha - 1.0, 2) + pow(beta - 1.0, 2));
      fit[j] = 1.0 - value[j];
      
    } else if (method == "NSE" || method == "Rsq") {
      double SS_res = 0;
      for(int i=0; i<n_obs; i++) {
        if(!NumericVector::is_na(obs[i]) && !NumericVector::is_na(cur_sim[i])) {
          SS_res += pow(cur_sim[i] - obs[i], 2);
        }
      }
      value[j] = 1.0 - (SS_res / SS_total);
      fit[j] = 1.0 - value[j];
      
    } else if (method == "Abias") {
      double diff_sum = 0;
      for(int i=0; i<n_obs; i++) {
        if(!NumericVector::is_na(obs[i]) && !NumericVector::is_na(cur_sim[i])) {
          diff_sum += (cur_sim[i] - obs[i]);
        }
      }
      value[j] = fabs((diff_sum / count_sim) / mean_obs);
      fit[j] = value[j];
    }
  }
  return List::create(Named("value") = value, Named("fit") = fit);
}