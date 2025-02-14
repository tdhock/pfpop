pfpop_map <- function(degrees_vec, penalty, weight_vec=rep(1,length(degrees_vec)), verbose_file=""){
  pfpop_map_interface(degrees_vec, penalty, weight_vec, verbose_file)
}
pfpop_list <- function(degrees_vec, penalty, weight_vec=rep(1,length(degrees_vec))){
  pfpop_list_interface(degrees_vec, penalty, weight_vec)
}
