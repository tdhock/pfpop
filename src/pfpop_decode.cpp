#include <stdio.h>

void pfpop_decode_internal
(const int *best_change_ptr,
 const double *best_param_ptr,
 const int N_data,
 const int N_segs,
 int *seg_start_ptr,
 int *seg_end_ptr,
 double *seg_param_ptr,
 int *seg_count
 ){
  int last_i = N_data-1;
  int seg_i = N_segs-1;
  if(N_segs==0)*seg_count = 0;
  while(0 <= last_i){
    int next_last = best_change_ptr[last_i];
    if(N_segs != 0){
      seg_start_ptr[seg_i] = next_last+1;
      seg_end_ptr[seg_i] = last_i;
      seg_param_ptr[seg_i] = best_param_ptr[last_i];
    }
    seg_i--;
    last_i = next_last;
    if(N_segs==0)(*seg_count)++;
  }
}  

void pfpop_decode
(const int *best_change_ptr,
 const double *best_param_ptr,
 const int N_data,
 const int N_segs,
 int *seg_start_ptr,
 int *seg_end_ptr,
 double *seg_param_ptr
 ){
  pfpop_decode_internal
    (best_change_ptr, best_param_ptr,
     N_data, N_segs,
     seg_start_ptr, seg_end_ptr, seg_param_ptr,
     0);
}

int pfpop_decode_size
(const int *best_change_ptr,
 const int N_data
 ){
  int N_segs;
  pfpop_decode_internal
    (best_change_ptr, 0,
     N_data, 0,
     0, 0, 0, &N_segs);
  return N_segs;
}
