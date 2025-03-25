#include <stdio.h>

int pfpop_decode
(const int *best_change_ptr,
 const double *best_cost_ptr,
 const double *best_param_ptr,
 const int N_data,
 int *seg_start_ptr,
 int *seg_end_ptr,
 double *seg_param_ptr,
 const int N_segs){
  int last_i = N_data-1;
  int seg_i = N_segs-1;
  int seg_count = 0;
  while(0 <= last_i){
    int next_last = best_change_ptr[last_i];
    if(N_segs != 0){
      seg_start_ptr[seg_i] = next_last+1;
      seg_end_ptr[seg_i] = last_i;
      seg_param_ptr[seg_i] = best_param_ptr[last_i];
    }
    seg_i--;
    last_i = next_last;
    seg_count++;
  }
  return seg_count;
}  
