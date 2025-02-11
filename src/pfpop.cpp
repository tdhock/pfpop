#include <iomanip> //for setprecision.
#include <fstream> //for ifstream etc.
#include <exception>//for std::exception
#include <stdexcept>//for std::invalid_argument
#include <R.h> // Rprintf

#include "pfpop.h"

#include <math.h>
#include <stdio.h>
#include <R.h>

#define PREV_NOT_SET (-1)
#define MAX_ANGLE 360

#define ABS(x) ((x)<0 ? -(x) : (x))

LinearCoefsForList::LinearCoefsForList
(double li, double co, double m, double M, int i, double prev){
  Linear = li;
  Constant = co;
  min_angle_param = m;
  max_angle_param = M;
  data_i = i;
  prev_angle_param = prev;
}

LinearCoefsForList::LinearCoefsForList
(double li, double co, double m, double M){
  Linear = li;
  Constant = co;
  min_angle_param = m;
  max_angle_param = M;
  data_i = PREV_NOT_SET;
  prev_angle_param = INFINITY;
}

LinearCoefsForList::LinearCoefsForList(){
}

double LinearCoefsForList::Loss(double angle_param){
  return Linear*angle_param + Constant;
}

L1LossListFun::L1LossListFun(){
  weight = 1;
}

void L1LossListFun::set_to_min_of_one
(L1LossListFun *input, int verbose){
  double best_loss = INFINITY,
    best_angle_param = INFINITY,
    prev_angle_param = INFINITY;
  int data_i;
  input->Minimize(&best_loss, &best_angle_param, &data_i);
  piece_list.clear();
  piece_list.emplace_front(0, best_loss, 0, MAX_ANGLE, PREV_NOT_SET, best_angle_param);
}

void L1LossListFun::push_sum_pieces
(L1LossListFun *fun1,
 L1LossListFun *fun2,
 L1LossList::iterator it1,
 L1LossList::iterator it2,
 int verbose){
  double min_angle_param =
    (it1->min_angle_param < it2->min_angle_param) ?
    it2->min_angle_param : it1->min_angle_param;
  double max_angle_param =
    (it1->max_angle_param < it2->max_angle_param) ?
    it1->max_angle_param : it2->max_angle_param;
  enlarge_last_or_emplace
    (it1->Linear+it2->Linear, it1->Constant+it2->Constant,
     min_angle_param, max_angle_param, it2->data_i, it2->prev_angle_param);
}

void L1LossListFun::push_min_pieces
(L1LossListFun *fun1,
 L1LossListFun *fun2,
 L1LossList::iterator it1,
 L1LossList::iterator it2,
 int verbose){
  double min_angle_param =
    (it1->min_angle_param < it2->min_angle_param) ?
    it2->min_angle_param : it1->min_angle_param;
  double max_angle_param =
    (it1->max_angle_param < it2->max_angle_param) ?
    it1->max_angle_param : it2->max_angle_param;
  double diff_Constant = it1->Constant-it2->Constant;
  double diff_Linear = it1->Linear-it2->Linear;
  double root_angle_param = -diff_Constant/diff_Linear;
  if(diff_Linear == 0){
    if(diff_Constant < 0){
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }else{
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }
  }else if(root_angle_param <= min_angle_param){
    if(diff_Linear < 0){
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }else{
      // f1-f2>0 => f1>f2
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }
  }else if(root_angle_param >= max_angle_param){
    if(diff_Linear < 0){
      // f1-f2>0 => f1>f2
      enlarge_last_or_emplace
	(it2, min_angle_param, max_angle_param);
    }else{
      // f1-f2<0 => f1<f2
      enlarge_last_or_emplace
	(it1, min_angle_param, max_angle_param);
    }
  }else{
    //it1 intersects it2 between min and max -> add two pieces.
    if(diff_Linear < 0){
      //f1-f2>0 => f1>f2 before.
      enlarge_last_or_emplace
	(it2, min_angle_param, root_angle_param);
      enlarge_last_or_emplace
	(it1, root_angle_param, max_angle_param);
    }else{
      enlarge_last_or_emplace
	(it1, min_angle_param, root_angle_param);
      enlarge_last_or_emplace
	(it2, root_angle_param, max_angle_param);
    }      
  }
}

void L1LossListFun::set_to_min_of_two
(L1LossListFun *fun1,
 L1LossListFun *fun2,
 int verbose){
  while_piece_pairs(fun1, fun2, &L1LossListFun::push_min_pieces, verbose);
}

void L1LossListFun::set_to_sum_of
(L1LossListFun *fun1,
 L1LossListFun *fun2,
 int verbose){
  while_piece_pairs(fun1, fun2, &L1LossListFun::push_sum_pieces, verbose);
}

void L1LossListFun::while_piece_pairs
(L1LossListFun *fun1,
 L1LossListFun *fun2,
 push_fun_ptr push_pieces,
 int verbose){
  L1LossList::iterator
    it1 = fun1->piece_list.begin(),
    it2 = fun2->piece_list.begin();
  piece_list.clear();
  while(it1 != fun1->piece_list.end() &&
	it2 != fun2->piece_list.end()){
    (this->*push_pieces)(fun1, fun2, it1, it2, verbose);
    double last_max_angle_param = piece_list.back().max_angle_param;
    if(it1->max_angle_param == last_max_angle_param){
      it1++;
    }
    if(it2->max_angle_param == last_max_angle_param){
      it2++;
    }
  }
}

void L1LossListFun::add(double Constant){
  L1LossList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Constant += Constant;
  }
}

void L1LossListFun::multiply(double x){
  L1LossList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->Linear *= x;
    it->Constant *= x;
  }
}

void L1LossListFun::set_prev_seg_end(int prev_seg_end){
  L1LossList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->data_i = prev_seg_end;
  }
}

void L1LossListFun::findMean
(double angle_param, int *seg_end, double *prev_angle_param){
  L1LossList::iterator it;
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    if(it->min_angle_param <= angle_param && angle_param <= it->max_angle_param){
      *seg_end = it->data_i;
      *prev_angle_param = it->prev_angle_param;
      return;
    }
  }
}

void L1LossListFun::print(){
  L1LossList::iterator it;
  Rprintf("%5s %5s %5s %5s %5s %s\n",
	  "Linear", "Constant",
	  "min_angle_param", "max_angle_param",
	  "data_i", "prev_angle_param");
  for(it=piece_list.begin(); it != piece_list.end(); it++){
    it->print();
  }
}

void LinearCoefsForList::print(){
  Rprintf("%.20e %.20e %15f %15f %15f %d\n",
	 Linear, Constant,
	 min_angle_param, max_angle_param,
	 prev_angle_param, data_i);
}

void L1LossListFun::emplace_piece
(double Linear, double Constant,
 double min_angle_param, double max_angle_param){
  piece_list.emplace_back
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param);
}

void L1LossListFun::emplace_piece
(double Linear, double Constant,
 double min_angle_param, double max_angle_param,
 int data_i, double prev_angle_param){
  piece_list.emplace_back
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param,
     data_i, prev_angle_param);
}

void L1LossListFun::enlarge_last_or_emplace
(double Linear, double Constant,
 double min_angle_param, double max_angle_param,
 int data_i, double prev_angle_param){
  L1LossList::iterator it=piece_list.end();
  if(it!=piece_list.begin()){
    it--;
    if(it->Linear == Linear && it->Constant == Constant &&
       it->data_i==data_i && it->prev_angle_param==prev_angle_param){
      it->max_angle_param = max_angle_param;
      return;
    }
  }
  emplace_piece
    (Linear*weight, Constant*weight,
     min_angle_param, max_angle_param,
     data_i, prev_angle_param);
}

void L1LossListFun::enlarge_last_or_emplace
(double Linear, double Constant,
 double min_angle_param, double max_angle_param){
  enlarge_last_or_emplace
    (Linear, Constant, min_angle_param, max_angle_param,
     PREV_NOT_SET, INFINITY);
}

void L1LossListFun::enlarge_last_or_emplace
(L1LossList::iterator it,
 double min_angle_param, double max_angle_param){
  enlarge_last_or_emplace
    (it->Linear, it->Constant,
     min_angle_param, max_angle_param,
     it->data_i, it->prev_angle_param);
}

void L1LossListFun::init
(double angle, double weight_){
  weight = weight_;
  piece_list.clear();
  if(angle == 0){
    emplace_piece(1, 0, 0, MAX_ANGLE/2);
    emplace_piece(-1, MAX_ANGLE, MAX_ANGLE/2, MAX_ANGLE);
  }else if(angle < MAX_ANGLE/2){
    emplace_piece(-1, angle, 0, angle);
    emplace_piece(1, -angle, angle, angle+MAX_ANGLE/2);
    emplace_piece(-1, (MAX_ANGLE+angle), angle+MAX_ANGLE/2, MAX_ANGLE);
  }else if(angle == MAX_ANGLE/2){
    emplace_piece(-1, MAX_ANGLE/2, 0, MAX_ANGLE/2);
    emplace_piece(1, -MAX_ANGLE/2, MAX_ANGLE/2, MAX_ANGLE);
  }else{
    emplace_piece(1, MAX_ANGLE-angle, 0, angle-MAX_ANGLE/2);
    emplace_piece(-1, angle, angle-MAX_ANGLE/2, angle);
    emplace_piece(1, -angle, angle, MAX_ANGLE);
  }
}

void L1LossListFun::Minimize
(double *best_loss,
 double *best_angle_param,
 int *data_i){
  *best_loss = INFINITY;
  for
    (L1LossList::iterator it = piece_list.begin();
     it != piece_list.end();
     it ++){
    double it_angle_param =
      (it->Linear < 0) ? it->max_angle_param : it->min_angle_param;
    double it_loss = it->Loss(it_angle_param);
    if(it_loss < *best_loss){
      *best_loss = it_loss;
      *best_angle_param = it_angle_param;
      *data_i = it->data_i;
    }
  }
}

int pfpop
(const double *degrees_ptr,
 const double penalty,
 const double *weight_ptr,
 const int N_data,
 int *best_change_ptr,
 double *best_cost_ptr,
 double *best_param_ptr,
 int *best_N_segs_ptr,
 int *num_pieces_ptr){
  if(penalty == INFINITY){
    //ok.
  }else if(!std::isfinite(penalty)){
    return ERROR_PENALTY_NOT_FINITE;
  }else if(penalty < 0){
    return ERROR_PENALTY_NEGATIVE;
  }
  L1LossListFun dist_fun_i, cost_up_to_i, cost_up_to_prev, cost_of_change, min_term;
  int verbose=0;
  double cum_weight_i = 0, cum_weight_prev_i = 0;
  double total_intervals = 0.0, max_intervals = 0.0;
  for(int data_i=0; data_i<N_data; data_i++){
    double angle = degrees_ptr[data_i];
    if(!std::isfinite(angle)){
      return ERROR_DATA_NOT_FINITE;
    }
    if(angle<0){
      return ERROR_DATA_NEGATIVE;
    }
    if(angle >= MAX_ANGLE){
      return ERROR_DATA_NOT_LESS_THAN_360;
    }
    double weight = weight_ptr[data_i];
    if(!std::isfinite(weight)){
      return ERROR_WEIGHT_NOT_FINITE;
    }
    if(weight <= 0){
      return ERROR_WEIGHT_NOT_POSITIVE;
    }
    cum_weight_i += weight;
    dist_fun_i.init(angle, weight);
    if(data_i==0){
      cost_up_to_i = dist_fun_i;
    }else{
      if(penalty == INFINITY){
	min_term = cost_up_to_prev;
      }else{
	cost_of_change.set_to_min_of_one(&cost_up_to_prev, verbose);
	// V_t(m) = (gamma_t + w_{1:t-1} * M_t(m))/w_{1:t}, where
	// M_t(m) = min{
	//   V_{t-1}(m),
	//   Vbar_{t-1} + penalty/w_{1:t-1}
	// in other words, we need to divide the penalty by the previous cumsum,
	// and add that to the min-less-ified function, before applying the min-env
	cost_of_change.set_prev_seg_end(data_i-1);
	cost_of_change.add(penalty/cum_weight_prev_i);
	if(penalty==0){
	  min_term = cost_of_change;
	}else{
	  min_term.set_to_min_of_two(&cost_of_change, &cost_up_to_prev, verbose);
	}
      }
      min_term.multiply(cum_weight_prev_i);
      cost_up_to_i.set_to_sum_of(&dist_fun_i, &min_term, verbose);
    }
    cost_up_to_i.multiply(1/cum_weight_i);
    cum_weight_prev_i = cum_weight_i;
    num_pieces_ptr[data_i] = cost_up_to_i.piece_list.size();
    cost_up_to_prev = cost_up_to_i;
    cost_up_to_i.Minimize
      (best_cost_ptr+data_i,
       best_param_ptr+data_i,
       best_change_ptr+data_i);
  }//while(can read line in text file)
  // Decoding the cost_model_vec, and writing to the output matrices.
  *best_N_segs_ptr = decode
    (best_change_ptr, best_cost_ptr, best_param_ptr, N_data,
     0, 0, 0, 0);
  return 0;
}

int decode
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
