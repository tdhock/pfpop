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
 const char *verbose_file,
 int *best_change_ptr,
 double *max_cost_ptr,
 double *max_param_ptr,
 double *max_Linear_ptr,
 double *max_Constant_ptr,
 double *min_cost_ptr,
 double *min_param_ptr,
 double *min_Linear_ptr,
 double *min_Constant_ptr,
 int *best_N_segs_ptr,
 int *map_size_ptr,
 int *list_size_ptr){
  std::ofstream verbose_fstream; // ofstream supports output only.
  bool verbose = strcmp(verbose_file, "") != 0;
  if(verbose){
    verbose_fstream.open(verbose_file);
    verbose_fstream << "data_i" << "\t" << "first_param" << "\t" << "opt_param" << "\t" << "last_param" << "\t" << "first_diff" << "\t" << "opt_diff" << "\t" << "last_diff" << "\t" << "Linear" << "\t" << "Constant" << "\n";
  }
  L1LossMapFun cost_model;
  cost_model.cost = 0;
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
    cost_model.write_min_or_max
      (data_i, -1,
       min_cost_ptr, min_param_ptr, min_Linear_ptr, min_Constant_ptr);
    cost_model.write_min_or_max
      (data_i, 1,
       max_cost_ptr, max_param_ptr, max_Linear_ptr, max_Constant_ptr);
    list_size_ptr[data_i] = cost_model.ptr_list.size();
    if(verbose){
      for
	(ClusterList::iterator cluster_it=cost_model.ptr_list.begin();
	 cluster_it != cost_model.ptr_list.end();
	 cluster_it++){
	verbose_fstream << data_i << "\t" <<
	  cost_model.get_param(cluster_it->first) << "\t" <<
	  cost_model.get_param(cluster_it->opt) << "\t" <<
	  cost_model.get_param(cluster_it->last) << "\t" <<
	  cost_model.get_Linear_diff(cluster_it->first) << "\t" <<
	  cost_model.get_Linear_diff(cluster_it->opt) << "\t" <<
	  cost_model.get_Linear_diff(cluster_it->last) << "\t" <<
	  cluster_it->opt.Linear << "\t" << cluster_it->opt.Constant << "\n";
      }
    }
  }
  best_N_segs_ptr[0]=0;//TODO.
  return 0;
}

double L1LossMapFun::get_param(Mapit &mit){
  return get_param(mit.it);
}

void L1LossMapFun::write_min_or_max
(int data_i, int sign,
 double *best, double *param, double *Linear, double *Constant
 ){
  best[data_i] = -INFINITY * sign;
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    double cost = get_cost_at_ptr(*it);
    if(sign*cost > sign*best[data_i]){
      best[data_i] = cost;
      if(param)param[data_i] = get_param_or_mid(*it);
      if(Linear)Linear[data_i] = it->opt.Linear;
      if(Constant)Constant[data_i] = it->opt.Constant;
    }
  }
}

int sgn(double x){
  if(x<0)return -1;
  if(x>0)return 1;
  return 0;
}

void L1LossMapFun::add_loss_for_data(double angle_, double weight_){
  angle = angle_;
  weight = weight_;
  step = 1;//add/update breaks
  pieces();
  step = 2;//update coefs
  all_pointers();
  step = 3;//move ptr if Linear_diff=0
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    move_right_if_zero(it->first);
    move_right_if_zero(it->opt);// points to piece on and after the breakpoint.
    move_left_if_zero(it->last);
  }
  step = 4;//delete break if Linear_diff=0
  pieces();
  step = 5;//combine adjacent pairs of pointers.
  // ClusterList::iterator it=ptr_list.begin(), next_it=ptr_list.begin()++;
  // while(next_it != ptr_list.end()){
  //   if(sgn(it->last.it->second.Linear_diff) == sgn(next_it->first.it->second.Linear_diff)){
  //     next_it->first.it = it->first.it;
  //     ptr_list.erase(it);
  //   }
  //   it = next_it;
  //   next_it++;
  // }
  //step = 6;//split pointers
  //all_pointers();
  step = 7;//move opt iterators.
  all_pointers();
}

void L1LossMapFun::all_pointers(){
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    cluster_it = it;
    pieces();
  }
}

double L1LossMapFun::max(){
  return min_or_max(1);
}

double L1LossMapFun::min(){
  return min_or_max(-1);
}

double L1LossMapFun::min_or_max(int sign){
  double best;
  write_min_or_max(0, -1, &best, 0, 0, 0);
  return best;
}

double L1LossMapFun::get_param_or_mid(const Cluster cl){
  if(cl.opt.Linear!=0)return(get_param(cl.opt.it));
  if(cl.opt.it == loss_map.end())return INFINITY;
  double last_param, first_param;
  if(cl.opt.it == loss_map.begin()){
    L1LossMap::iterator last_it = loss_map.end();
    last_it--;
    last_param = get_param(last_it);
    first_param = get_param(cl.opt.it) + MAX_ANGLE;
  }else{
    L1LossMap::iterator prev_it = cl.opt.it;
    prev_it--;
    last_param = get_param(cl.opt.it);
    first_param = get_param(prev_it);
  }
  return (last_param+first_param)/2;
}

double L1LossMapFun::get_cost_at_ptr(const Cluster cl){
  return get_param(cl.opt.it)*cl.opt.Linear+cl.opt.Constant;
}

bool between(double first, double param, double last){
  double min, max;
  if(first < last){
    min = first;
    max = last;
  }else{
    min = last;
    max = first;
  }
  bool in_min_max = min < param && param < max;
  return (first < last) ? in_min_max : !in_min_max;
}

void L1LossMapFun::piece
(double Linear_, double Constant_,
 double min_param_, double max_param_){
  Linear=Linear_*weight;
  Constant=Constant_*weight;
  min_param=min_param_;
  max_param=max_param_;
  double diff_Linear_at_min = (min_param==0 && max_param != MAX_ANGLE/2) ? 0 : 2*Linear;
  if(step==1 && diff_Linear_at_min){//insert or update bkpt in map
    std::pair<L1LossMap::iterator, bool> result;
    Breakpoint new_break(0, -2);
    result = loss_map.insert(std::pair<double,Breakpoint>(min_param, new_break));
    L1LossMap::iterator insert_it = result.first;
    if(result.second){
      // new breakpoint inserted, so we need to copy the adjacent data_i.
      L1LossMap::iterator adj_it = insert_it;
      adj_it++;
      set_data_i(insert_it, get_data_i(adj_it));
    }
    insert_it->second.Linear_diff += diff_Linear_at_min;
    if(ptr_list.size()<2){
      Cluster new_cl;
      new_cl.first.it = new_cl.opt.it = new_cl.last.it = insert_it;
      new_cl.sign = sgn(insert_it->second.Linear_diff);
      ptr_list.push_back(new_cl);
    }
  }
  if(step==2 && min_param <= get_param(cluster_it->opt) && get_param(cluster_it->opt) < max_param){//update coefs
    cluster_it->opt.Linear += Linear;
    cluster_it->opt.Constant += Constant;
  }
  if(step==4 && diff_Linear_at_min){//delete break if diff_linear=0.
    L1LossMap::iterator it = loss_map.find(min_param);
    if(it->second.Linear_diff==0){
      loss_map.erase(it);
    }
    // TODO adjust pointers if sign changed on edge!!
  }
  if
    (step==6 &&//split if necessary
     diff_Linear_at_min &&
     between(cluster_it->first.it->first, min_param, cluster_it->last.it->first)
     ){
    L1LossMap::iterator it = loss_map.find(min_param);
    if(it!=loss_map.end() && cluster_it->sign != sgn(it->second.Linear_diff)){
      Cluster new_cl;
      new_cl.sign = -cluster_it->sign;
      new_cl.opt = cluster_it->opt;
      move_it_fun_ptr move;
      bool moving_left = between
	(cluster_it->first.it->first,
	 min_param,
	 cluster_it->opt.it->first);
      ClusterList::iterator insert_it = cluster_it;
      if(moving_left){
	move = &L1LossMapFun::move_left;
      }else{
	move = &L1LossMapFun::move_right;
	insert_it++;
      }
      while(new_cl.opt.it != it){
	(this->*move)(new_cl.opt);
      }
      new_cl.first.it = cluster_it->opt.it;
      new_cl.last.it = cluster_it->opt.it;
      new_cl.optimize();
      ClusterList::iterator new_it = ptr_list.insert(cluster_it, new_cl);
      if(moving_left){
	new_cl.first.it = cluster_it->first.it;
	cluster_it->first.it = new_cl.opt.it;
	cluster_it->first.it++;
	insert_it = new_it;
      }else{
	new_cl.last.it = cluster_it->first.it;
	cluster_it->last.it = new_cl.opt.it;
	cluster_it->last.it--;
      }
      new_cl.optimize();
      new_cl.sign = -cluster_it->sign;
      ptr_list.insert(insert_it, new_cl);
    }
  }
  if(step==7){//move opt it
  }
}

void Cluster::optimize(){
  //TODO
}

void L1LossMapFun::move_right_if_zero(Mapit &mit){
  move_if_zero(&L1LossMapFun::move_right, mit);
}
void L1LossMapFun::move_left_if_zero(Mapit &mit){
  move_if_zero(&L1LossMapFun::move_left, mit);
}
void L1LossMapFun::move_if_zero(move_it_fun_ptr move, Mapit &mit){
  if(mit.it != loss_map.end() && mit.it->second.Linear_diff==0){
    (this->*move)(mit);
    if(mit.it->second.Linear_diff==0){
      // if after moving we are still at a zero, they all must be
      // zeros, so move to the end.
      mit.it = loss_map.end();
    }
  }
}

void L1LossMapFun::move_left(Mapit &mit){
  update_coefs(mit, -1, get_param(mit), get_param(mit));
  if(mit.it == loss_map.begin()){
    update_coefs(mit, -1, 0, MAX_ANGLE);
    mit.it = loss_map.end();
  }
  mit.it--;
}

void L1LossMapFun::move_right(Mapit &mit){
  mit.it++;
  if(mit.it == loss_map.end()){
    update_coefs(mit, 1, MAX_ANGLE, 0);
    mit.it = loss_map.begin();
  }
  update_coefs(mit, 1, get_param(mit), get_param(mit));
}

void L1LossMapFun::update_coefs
(Mapit &mit, int sign, double min_param, double max_param){
  mit.update_coefs(sign, min_param, max_param, get_Linear_diff(mit.it));
}

void Mapit::update_coefs
(int sign, double param_before, double param_after, double ldiff){
  // do nothing.
}

void Coefs::update_coefs
(int sign, double param_before, double param_after, double ldiff){
  double cost = Constant + Linear*param_before;
  Linear += sign*ldiff;
  Constant = cost-Linear*param_after;
}

double L1LossMapFun::get_param(L1LossMap::iterator it){
  if(it == loss_map.end())return MAX_ANGLE;
  return it->first;
}

void L1LossMapFun::pieces(){
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
  // TODO initial Clusters?
}

Breakpoint::Breakpoint(double l,int d){
  Linear_diff=l;
  data_i=d;
}

Breakpoint::Breakpoint(){
  Linear_diff=0;
  data_i=-1;
}

int L1LossMapFun::get_data_i(L1LossMap::iterator it){
  if(it == loss_map.end())return -9;
  return it->second.data_i;
}

double L1LossMapFun::get_Linear_diff(Mapit &mit){
  return get_Linear_diff(mit.it);
}

double L1LossMapFun::get_Linear_diff(L1LossMap::iterator it){
  if(it == loss_map.end())return INFINITY;
  return it->second.Linear_diff;
}

void L1LossMapFun::set_data_i(L1LossMap::iterator it, int data_i){
  if(it == loss_map.end())return;
  it->second.data_i = data_i;
}

void L1LossMapFun::min_with_constant(double constant){
}

L1LossMapFun::L1LossMapFun(){
  delete_breaks(0);
}

Cluster::Cluster(){
}

void Cluster::init(L1LossMap::iterator it_, double Constant_){
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
