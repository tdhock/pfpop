#include "pfpop.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List pfpop_interface
(const Rcpp::NumericVector degrees_vec,
 const double penalty,
 const Rcpp::NumericVector weight_vec){
  int N_data = degrees_vec.length();
  if(N_data < 1){
    Rcpp::stop("degrees_vec length must be one or more");
  }
  if(N_data != weight_vec.length()){
    Rcpp::stop("degrees_vec and weight_vec lengths must be equal");
  }
  Rcpp::IntegerVector best_change_vec(N_data);
  Rcpp::NumericVector best_cost_vec(N_data);
  Rcpp::NumericVector best_param_vec(N_data);
  Rcpp::IntegerVector best_N_segs(1);
  Rcpp::IntegerVector num_pieces_vec(N_data);
  int status = pfpop_list
    (&degrees_vec[0],
     penalty,
     &weight_vec[0],
     N_data,
     &best_change_vec[0],
     &best_cost_vec[0],
     &best_param_vec[0],
     &best_N_segs[0],
     &num_pieces_vec[0]);
  if(status==ERROR_PENALTY_NEGATIVE || status==ERROR_PENALTY_NOT_FINITE){
    Rcpp::stop("penalty=%f must be non-negative", penalty);
  }
  if(status==ERROR_DATA_NEGATIVE || status==ERROR_DATA_NOT_LESS_THAN_360 || status==ERROR_DATA_NOT_FINITE){
    Rcpp::stop("data values must be in [0,360)");
  }
  if(status==ERROR_WEIGHT_NOT_POSITIVE || status==ERROR_WEIGHT_NOT_FINITE){
    Rcpp::stop("weight values must be positive");
  }
  if(status != 0){
    Rcpp::stop("error code %d", status);
  }
  int N_segs = best_N_segs[0];
  Rcpp::IntegerVector seg_start_vec(N_segs);
  Rcpp::IntegerVector seg_end_vec(N_segs);
  Rcpp::NumericVector seg_param_vec(N_segs);
  int ignored = decode
    (&best_change_vec[0],
     &best_cost_vec[0],
     &best_param_vec[0],
     N_data,
     &seg_start_vec[0],
     &seg_end_vec[0],
     &seg_param_vec[0],
     N_segs);
  return Rcpp::List::create
    (Rcpp::Named
     ("segments", Rcpp::DataFrame::create
      (Rcpp::Named("start", seg_start_vec),
       Rcpp::Named("end", seg_end_vec),
       Rcpp::Named("param", seg_param_vec))),
     Rcpp::Named
     ("iterations", Rcpp::DataFrame::create
      (Rcpp::Named("change", best_change_vec),
       Rcpp::Named("cost", best_cost_vec),
       Rcpp::Named("param", best_param_vec),
       Rcpp::Named("pieces", num_pieces_vec))));
}
