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

int pfpop_map
(const double *degrees_ptr,
 const double penalty,
 const double *weight_ptr,
 const int N_data,
 int *best_change_ptr,
 double *best_cost_ptr,
 double *best_param_ptr,
 int *best_N_segs_ptr,
 int *map_size_ptr,
 int *pointer_moves_ptr){
  L1LossMapFun cost_model;
  double cum_weight_i = 0, cum_weight_prev_i = 0;
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
    if(data_i != 0){
      if(penalty <= cost_model.min()){
	cost_model.delete_breaks(penalty);
      }else if(cost_model.max() <= penalty){
	// no update to cost model.
      }else{
	//cost_model.min_with_constant(penalty);
	return ERROR_NOT_IMPLEMENTED;
      }
    }
    // TODO to compute the mean cost instead of the total cost, we
    // divide the penalty by the previous cumsum, and add that to the
    // min-ified constant, before applying the min with constant.
    cost_model.add_loss_for_data(angle, weight);
    pointer_moves_ptr[data_i] = cost_model.move_both_pointers();
  }
  return 0;
}

int L1LossMapFun::move_both_pointers(){
  int moves = 0;
  moves += move_one_pointer(max_ptr, 1);
  moves += move_one_pointer(min_ptr, -1);
  return moves;
}

int L1LossMapFun::move_one_pointer(Pointer &ptr, int sign){
  if(ptr.Linear * sign < 0){
    move_left(ptr);
    return 1;
  }
  if(next_Linear(ptr) * sign >= 0){
    move_right(ptr);
    return 1;
  }
  return 0;
}

double L1LossMapFun::next_Linear(Pointer &ptr){
  return ptr.Linear + it->second->Linear_diff;
}

void L1LossMapFun::move_left(Pointer &ptr){
  if(ptr.it == loss_map.begin()){
    ptr.it = loss_map.end();
    Breakpoint *bptr = get_break_ptr(ptr.it);
    if(bptr->Linear_diff == 0){
      ptr.it--;
    }
  }
  TODO;
  update_coefs(ptr, -1);
}

void L1LossMapFun::move_right(Pointer &ptr){
  TODO;
  update_coefs(1);
}

double L1LossMapFun::min(){
  return get_cost_at_pointer(min_ptr);
}

double L1LossMapFun::max(){
  return get_cost_at_pointer(max_ptr);
}

double L1LossMapFun::get_cost_at_pointer(Pointer& pointer){
  return get_param(pointer.it)*Linear+Constant;
}

void L1LossMapFun::piece
(double Linear_, double Constant_,
 double min_param_, double max_param_){
  Linear=Linear_*weight;
  Constant=Constant_*weight;
  min_param=min_param_;
  max_param=max_param_;
  maybe_add(&min_ptr);
  maybe_add(&max_ptr);
  Breakpoint* break_ptr;
  if(max_param==MAX_ANGLE){
    break_ptr = &last_break;
  }else{
    //https://cplusplus.com/reference/map/map/emplace/ says If the
    //function successfully inserts the element (because no equivalent
    //element existed already in the map), the function returns a pair
    //of an iterator to the newly inserted element and a value of
    //true. Otherwise, it returns an iterator to the equivalent
    //element within the container and a value of false.
    std::pair<L1LossMap::iterator, bool> result;
    double diff_Linear = (max_param==MAX_ANGLE && min_param!=MAX_ANGLE/2) ? 0 : -2*Linear;
    Breakpoint new_break(0, -2);
    result = loss_map.insert(std::pair<double,Breakpoint>(max_param, new_break));
    //result = loss_map.emplace(max_param, INFINITY, -1);
    result.first->second.Linear_diff += diff_Linear;
    if(result.second){
      // new breakpoint inserted, so we need to copy the adjacent data_i.
      L1LossMap::iterator adj_it = result.first;
      adj_it++;
      set_data_i(result.first, get_data_i(adj_it));
    }
  }
  //TODO add coefs;
  //TODO remove break if necessary;
}

double L1LossMapFun::get_param(L1LossMap::iterator it){
  if(it == loss_map.end())return MAX_ANGLE;
  return it->first;
}

void L1LossMapFun::maybe_add(Pointer* ptr){
  double param = get_param(ptr->it);
  if(min_param < param & param <= max_param){
    ptr->Constant += Constant;
    ptr->Linear += Linear;
  }
}

void L1LossMapFun::add_loss_for_data(double angle, double weight_){
  weight = weight_;
  if(angle == 0){
    piece(1, 0, 0, MAX_ANGLE/2);
    piece(-1, MAX_ANGLE, MAX_ANGLE/2, MAX_ANGLE);
  }else if(angle < MAX_ANGLE/2){
    piece(-1, angle, 0, angle);
    piece(1, -angle, angle, angle+MAX_ANGLE/2);
    piece(-1, (MAX_ANGLE+angle), angle+MAX_ANGLE/2, MAX_ANGLE);
  }else if(angle == MAX_ANGLE/2){
    piece(-1, MAX_ANGLE/2, 0, MAX_ANGLE/2);
    piece(1, -MAX_ANGLE/2, MAX_ANGLE/2, MAX_ANGLE);
  }else{
    piece(1, MAX_ANGLE-angle, 0, angle-MAX_ANGLE/2);
    piece(-1, angle, angle-MAX_ANGLE/2, angle);
    piece(1, -angle, angle, MAX_ANGLE);
  }
}

void L1LossMapFun::delete_breaks(double new_cost){
  loss_map.clear();
  min_ptr.init(loss_map.begin(),new_cost);
  max_ptr.init(loss_map.begin(),new_cost);
}

Breakpoint::Breakpoint(double l,int d){
  Linear_diff=l;
  data_i=d;
}

Breakpoint::Breakpoint(){
  Linear_diff=0;
  data_i=-1;
}

Breakpoint* L1LossMapFun::get_break_ptr(L1LossMap::iterator it){
  if(it == loss_map.end())return &last_break;
  return &(it->second);
}

int L1LossMapFun::get_data_i(L1LossMap::iterator it){
  return get_break_ptr(it)->data_i;
}

double L1LossMapFun::get_Linear_diff(L1LossMap::iterator it){
  return get_break_ptr(it)->Linear_diff;
}

void L1LossMapFun::set_data_i(L1LossMap::iterator it, int data_i){
  get_break_ptr(it)->data_i = data_i;
}

void L1LossMapFun::set_Linear_diff(L1LossMap::iterator it, double Linear_diff){
  get_break_ptr(it)->Linear_diff = Linear_diff;
}

void L1LossMapFun::min_with_constant(double constant){
}

L1LossMapFun::L1LossMapFun(){
  delete_breaks(0);
}

Pointer::Pointer(){
  Constant=0;
  Linear=0;
}

void Pointer::init(L1LossMap::iterator it_, double Constant_){
  Constant=Constant_;
  Linear=0;
  it = it_;
}

int pfpop_list
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
