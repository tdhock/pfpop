#include "pfpop_map_l1.h"
#include "pfpop_list_l1.h"
#include "pfpop_decode.h"
#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List pfpop_map_l1_interface
(const Rcpp::NumericVector data_vec,
 const double penalty,
 const Rcpp::NumericVector weight_vec,
 const std::string verbose_file){
  int N_data = data_vec.length();
  if(N_data < 1){
    Rcpp::stop("data_vec length must be one or more");
  }
  if(N_data != weight_vec.length()){
    Rcpp::stop("data_vec and weight_vec lengths must be equal");
  }
  Rcpp::NumericVector cost_vec(N_data);
  Rcpp::NumericVector param_vec(N_data);
  Rcpp::IntegerVector change_vec(N_data);
  Rcpp::IntegerVector map_size_vec(N_data);
  Rcpp::IntegerVector list_size_vec(N_data);
  Rcpp::IntegerVector num_moves_vec(N_data);
  int status = pfpop_map_l1
    (&data_vec[0],
     penalty,
     &weight_vec[0],
     N_data,
     verbose_file.c_str(),
     //Inputs above, outputs below.
     &cost_vec[0],
     &param_vec[0],
     &change_vec[0],
     &map_size_vec[0],
     &list_size_vec[0],
     &num_moves_vec[0]);
  if(status==pfpop_map_ERROR_PENALTY_NEGATIVE || status==pfpop_map_ERROR_PENALTY_NOT_FINITE){
    Rcpp::stop("penalty=%f must be non-negative", penalty);
  }
  if(status==pfpop_map_ERROR_WEIGHT_NOT_POSITIVE || status==pfpop_map_ERROR_WEIGHT_NOT_FINITE){
    Rcpp::stop("weight values must be positive");
  }
  if(status != 0){
    Rcpp::stop("error code %d", status);
  }
  // Decoding the cost_model_vec, and writing to the output matrices.
  int N_segs = pfpop_decode_size(&change_vec[0], N_data);
  Rcpp::IntegerVector seg_start_vec(N_segs);
  Rcpp::IntegerVector seg_end_vec(N_segs);
  Rcpp::NumericVector seg_param_vec(N_segs);
  pfpop_decode
    (&change_vec[0],
     &param_vec[0],
     N_data,
     N_segs,
     &seg_start_vec[0],
     &seg_end_vec[0],
     &seg_param_vec[0]);
  return Rcpp::List::create
    (Rcpp::Named
     ("segments", Rcpp::DataFrame::create
      (Rcpp::Named("start", seg_start_vec),
       Rcpp::Named("end", seg_end_vec),
       Rcpp::Named("param", seg_param_vec))),
     Rcpp::Named
     ("iterations", Rcpp::DataFrame::create
      (Rcpp::Named("cost", cost_vec),
       Rcpp::Named("param", param_vec),
       Rcpp::Named("change", change_vec),
       Rcpp::Named("map_size", map_size_vec),
       Rcpp::Named("list_size", list_size_vec),
       Rcpp::Named("num_moves", num_moves_vec))));
}

// [[Rcpp::export]]
Rcpp::List pfpop_list_l1_interface
(const Rcpp::NumericVector data_vec,
 const double penalty,
 const Rcpp::NumericVector weight_vec,
 const std::string verbose_file){
  int N_data = data_vec.length();
  if(N_data < 1){
    Rcpp::stop("data_vec length must be one or more");
  }
  if(N_data != weight_vec.length()){
    Rcpp::stop("data_vec and weight_vec lengths must be equal");
  }
  Rcpp::IntegerVector best_change_vec(N_data);
  Rcpp::NumericVector best_cost_vec(N_data);
  Rcpp::NumericVector best_param_vec(N_data);
  Rcpp::IntegerVector num_pieces_vec(N_data);
  int status = pfpop_list_l1
    (&data_vec[0],
     penalty,
     &weight_vec[0],
     N_data,
     verbose_file.c_str(),
     //inputs above, outputs below.
     &best_change_vec[0],
     &best_cost_vec[0],
     &best_param_vec[0],
     &num_pieces_vec[0]);
  if(status==pfpop_list_ERROR_PENALTY_NEGATIVE || status==pfpop_list_ERROR_PENALTY_NOT_FINITE){
    Rcpp::stop("penalty=%f must be non-negative", penalty);
  }
  if(status==pfpop_list_ERROR_DATA_NEGATIVE || status==pfpop_list_ERROR_DATA_NOT_LESS_THAN_360 || status==pfpop_list_ERROR_DATA_NOT_FINITE){
    Rcpp::stop("data values must be in [0,360)");
  }
  if(status==pfpop_list_ERROR_WEIGHT_NOT_POSITIVE || status==pfpop_list_ERROR_WEIGHT_NOT_FINITE){
    Rcpp::stop("weight values must be positive");
  }
  if(status != 0){
    Rcpp::stop("error code %d", status);
  }
  // Decoding the cost_model_vec, and writing to the output matrices.
  int N_segs = pfpop_decode_size(&best_change_vec[0], N_data);
  Rcpp::IntegerVector seg_start_vec(N_segs);
  Rcpp::IntegerVector seg_end_vec(N_segs);
  Rcpp::NumericVector seg_param_vec(N_segs);
  pfpop_decode
    (&best_change_vec[0],
     &best_param_vec[0],
     N_data,
     N_segs,
     &seg_start_vec[0],
     &seg_end_vec[0],
     &seg_param_vec[0]);
  return Rcpp::List::create
    (Rcpp::Named
     ("segments", Rcpp::DataFrame::create
      (Rcpp::Named("start", seg_start_vec),
       Rcpp::Named("end", seg_end_vec),
       Rcpp::Named("param", seg_param_vec))),
     Rcpp::Named
     ("iterations", Rcpp::DataFrame::create
      (Rcpp::Named("cost", best_cost_vec),
       Rcpp::Named("param", best_param_vec),
       Rcpp::Named("change", best_change_vec),
       Rcpp::Named("num_pieces", num_pieces_vec))));
}
