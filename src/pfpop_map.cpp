#include <iomanip> //for setprecision.
#include <fstream> //for ifstream etc.
#include <exception>//for std::exception
#include <stdexcept>//for std::invalid_argument
#include <R.h> // Rprintf
#include <math.h>
#include <stdio.h>

#include "pfpop_map.h"

#define PREV_NOT_SET (-1)
#define MAX_ANGLE 360

class ClusterWriter {
public:
  std::ofstream verbose_fstream, breaks_fstream;
  L1LossMapFun *cost_model_ptr;
  void write(int data_i, int step_i){
    for
      (L1LossMap::iterator it=cost_model_ptr->loss_map.begin();
       it != cost_model_ptr->loss_map.end();
       it++){
      breaks_fstream << data_i << "\t" << step_i << "\t" << it->first << "\t" << it->second << "\n";
    }
    for
      (ClusterList::iterator cluster_it=cost_model_ptr->ptr_list.begin();
       cluster_it != cost_model_ptr->ptr_list.end();
       cluster_it++){
      verbose_fstream << data_i << "\t" <<
	step_i << "\t" << 
	cost_model_ptr->get_param(cluster_it->first) << "\t" <<
	cost_model_ptr->get_param(cluster_it->opt) << "\t" <<
	cost_model_ptr->get_param(cluster_it->last) << "\t" <<
	cluster_it->first.Linear << "\t" <<
	cluster_it->opt.Linear << "\t" <<
	cluster_it->last.Linear << "\t" <<
	cluster_it->first.Constant << "\t" <<
	cluster_it->opt.Constant << "\t" <<
	cluster_it->last.Constant << "\t" <<
	cluster_it->sign << "\n";
    }
  }
};

int pfpop_map
(const double *degrees_ptr,
 const double penalty,
 const double *weight_ptr,
 const int N_data,
 const char *verbose_file,
 double *max_cost_ptr,
 double *max_param_ptr,
 double *max_Linear_ptr,
 double *max_Constant_ptr,
 int *argmax_ptr,
 double *min_cost_ptr,
 double *min_param_ptr,
 double *min_Linear_ptr,
 double *min_Constant_ptr,
 int *argmin_ptr,
 int *map_size_ptr,
 int *list_size_ptr,
 int *num_moves_ptr){
  ClusterWriter vwriter;
  bool verbose = strcmp(verbose_file, "") != 0;
  std::string breaks_file = verbose_file;
  breaks_file += "_breaks";
  L1LossMapFun cost_model;
  if(verbose){
    vwriter.breaks_fstream.open(breaks_file);
    vwriter.breaks_fstream << "data_i" << "\t" << "step_i" << "\t" << "param" << "\t" << "Linear_diff" << "\n";
    vwriter.verbose_fstream.open(verbose_file);
    vwriter.verbose_fstream << "data_i" << "\t" << "step_i" << "\t" << "first_param" << "\t" << "opt_param" << "\t" << "last_param" << "\t" << "first_Linear" << "\t" << "opt_Linear" << "\t" << "last_Linear" << "\t" << "first_Constant" << "\t" << "opt_Constant" << "\t" << "last_Constant" << "\t" << "sign" << "\n";
    vwriter.cost_model_ptr = &cost_model;
  }
  cost_model.cost = 0;
  double cum_weight_i = 0, cum_weight_prev_i = 0;
  for(int data_i=0; data_i<N_data; data_i++){
    double angle = degrees_ptr[data_i];
    if(!std::isfinite(angle)){
      return pfpop_map_ERROR_DATA_NOT_FINITE;
    }
    if(angle<0){
      return pfpop_map_ERROR_DATA_NEGATIVE;
    }
    if(angle >= MAX_ANGLE){
      return pfpop_map_ERROR_DATA_NOT_LESS_THAN_360;
    }
    double weight = weight_ptr[data_i];
    if(!std::isfinite(weight)){
      return pfpop_map_ERROR_WEIGHT_NOT_FINITE;
    }
    if(weight <= 0){
      return pfpop_map_ERROR_WEIGHT_NOT_POSITIVE;
    }
    cum_weight_i += weight;
    cost_model.moves = 0;
    cost_model.data_i = data_i;
    printf("------>data_i=%d\n", data_i);
    if(data_i != 0){
      cost_model.min_with_constant(min_cost_ptr[data_i-1]+penalty);
    }
    if(verbose)vwriter.write(data_i, 0);
    // TODO to compute the mean cost instead of the total cost, we
    // divide the penalty by the previous cumsum, and add that to the
    // min-ified constant, before applying the min with constant.
    cost_model.add_loss_for_data(angle, weight);
    cost_model.write_min_or_max
      (data_i, -1,
       min_cost_ptr, min_param_ptr, min_Linear_ptr, min_Constant_ptr, argmin_ptr);
    cost_model.write_min_or_max
      (data_i, 1,
       max_cost_ptr, max_param_ptr, max_Linear_ptr, max_Constant_ptr, argmax_ptr);
    map_size_ptr[data_i] = cost_model.loss_map.size();
    list_size_ptr[data_i] = cost_model.ptr_list.size();
    num_moves_ptr[data_i] = cost_model.moves;
    if(verbose)vwriter.write(data_i, 1);
  }
  return 0;
}

double L1LossMapFun::get_param(Coefs &mit){
  return get_param(mit.it);
}

void L1LossMapFun::write_min_or_max
(int data_i, int sign,
 double *best, double *param, double *Linear, double *Constant, int *arg
 ){
  best[data_i] = -INFINITY * sign;
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    double cost = get_cost_at_coefs(it->opt);
    if(sign*cost > sign*best[data_i]){
      best[data_i] = cost;
      if(param)param[data_i] = get_param_or_mid(*it);
      if(Linear)Linear[data_i] = it->opt.Linear;
      if(Constant)Constant[data_i] = it->opt.Constant;
      if(arg)arg[data_i] = it->data_i;
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
    double first_cost = get_cost_at_coefs(it->first);
    double opt_cost = get_cost_at_coefs(it->opt);
    double last_cost = get_cost_at_coefs(it->last);
    //printf("add first=%f opt=%f last=%f\n", first_cost, opt_cost, last_cost);
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
  step = 6;//split pointers
  all_pointers();
  //step = 7;//move opt iterators.
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    double first_cost, opt_cost, last_cost;
    first_cost = get_cost_at_coefs(it->first);
    opt_cost = get_cost_at_coefs(it->opt);
    last_cost = get_cost_at_coefs(it->last);
    //printf("mid first=%f opt=%f last=%f\n", first_cost, opt_cost, last_cost);
    move_left(it->first);
    if(sgn(get_Linear_diff(it->first))!=it->sign){
      move_right(it->first);
    }
    //printf("%f before moving last sign=%d\n", get_Linear_diff(it->last), it->sign);
    move_right(it->last);
    //printf("%f after move right sign=%d\n", get_Linear_diff(it->last), it->sign);
    if(sgn(get_Linear_diff(it->last))!=it->sign){
      //printf("moving last left\n");
      move_left(it->last);
    }
    //printf("%f after moving last sign=%d\n", get_Linear_diff(it->last), it->sign);
    first_cost = get_cost_at_coefs(it->first);
    opt_cost = get_cost_at_coefs(it->opt);
    last_cost = get_cost_at_coefs(it->last);
    //printf("before move_to_opt first=%f opt=%f last=%f\n", first_cost, opt_cost, last_cost);
    move_to_opt(it);
  }
}

void L1LossMapFun::move_to_opt(ClusterList::iterator &it){
  while(prev_Linear(it->opt) * it->sign >= 0 && it->opt.it != it->first.it){
    //printf("Linear=%f prev=%f sign=%d move opt left\n", it->opt.Linear, prev_Linear(it->opt), it->sign);
    move_left(it->opt);
  }
  while(it->opt.Linear * it->sign < 0 && it->opt.it != it->last.it){
    //printf("move opt right\n");
    move_right(it->opt);
  }
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

double L1LossMapFun::get_param_or_mid(const Cluster cl){
  return(get_param(cl.opt.it));
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

double L1LossMapFun::get_cost_at_coefs(const Coefs coefs){
  return get_param(coefs.it)*coefs.Linear+coefs.Constant;
}

bool cost_between(double x, double m, double y){
  return (x<m && m<y) || (y<m && m<x);
}

bool param_between(double first, double param, double last){
  if(first==last)return false;
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
    result = loss_map.insert(std::pair<double,double>(min_param, 0));
    L1LossMap::iterator insert_it = result.first;
    insert_it->second += diff_Linear_at_min;
    if(ptr_list.size()<2){
      Cluster new_cl;
      new_cl.opt.Linear = 0;
      new_cl.opt.Constant = cost;
      new_cl.opt.it = insert_it;
      new_cl.first = new_cl.last = new_cl.opt;
      new_cl.sign = sgn(insert_it->second);
      new_cl.data_i = data_i-1;//?
      ptr_list.push_back(new_cl);
    }
  }
  if(step==2){
    update_coefs(cluster_it->first);
    update_coefs(cluster_it->opt);
    update_coefs(cluster_it->last);
  }
  if(step==4 && diff_Linear_at_min){//delete break if diff_linear=0.
    L1LossMap::iterator it = loss_map.find(min_param);
    if(it->second==0){
      loss_map.erase(it);
    }
    // TODO adjust pointers if sign changed on edge!!
  }
  //split if necessary
  if(step==6 && param_between
     (get_param(cluster_it->first),
      min_param,
      get_param(cluster_it->last))){
    L1LossMap::iterator it = loss_map.find(min_param);
    if(it!=loss_map.end() && // TODO not necessary?
       cluster_it->sign != sgn(it->second)){
      Cluster new_cl = *cluster_it;
      new_cl.sign = -cluster_it->sign;
      move_it_fun_ptr move;
      bool moving_left = param_between
	(get_param(cluster_it->first),
	 min_param,
	 get_param(cluster_it->opt));
      ClusterList::iterator insert_it = cluster_it;
      Coefs orig_end;
      if(moving_left){
	move = &L1LossMapFun::move_left;
	orig_end = cluster_it->first;
      }else{
	move = &L1LossMapFun::move_right;
	orig_end = cluster_it->last;
	insert_it++;
      }
      //printf("before while\n");
      while(new_cl.opt.it != it){
	if(moving_left){
	  cluster_it->first = new_cl.opt;
	}else{
	  cluster_it->last = new_cl.opt;
	}
	(this->*move)(new_cl.opt);
      }
      //printf("after while\n");
      new_cl.first = new_cl.opt;
      new_cl.last = new_cl.opt;
      //new_cl.optimize();
      ptr_list.insert(insert_it, new_cl);
      //printf("first new first=%f opt=%f last=%f\n", get_param(new_cl.first), get_param(new_cl.opt), get_param(new_cl.last));
      //printf("made it past insert\n");
      (this->*move)(new_cl.opt);//keep going one move past the new bkpt.
      if(moving_left){
	new_cl.first = orig_end;
	new_cl.last = new_cl.opt;
      }else{
	new_cl.first = new_cl.opt;
	new_cl.last = orig_end;
      }
      move_to_opt(cluster_it);
      //printf("cluster_it first=%f opt=%f last=%f\n", get_param(cluster_it->first), get_param(cluster_it->opt), get_param(cluster_it->last));
      new_cl.sign = cluster_it->sign;
      ClusterList::iterator new_it = ptr_list.insert(insert_it, new_cl);
      move_to_opt(new_it);
      //printf("second new first=%f opt=%f last=%f\n", get_param(new_cl.first), get_param(new_cl.opt), get_param(new_cl.last));
    }
  }
}

void L1LossMapFun::update_coefs(Coefs &coefs){
  if(min_param <= get_param(coefs) && get_param(coefs) < max_param){
    printf("updating param=%f Linear=%f+%f Constant=%f+%f\n", get_param(coefs), coefs.Linear, Linear, coefs.Constant, Constant);
    coefs.Linear += Linear;
    coefs.Constant += Constant;
  }
}

double L1LossMapFun::prev_Linear(Coefs &mit){
  return mit.Linear - get_Linear_diff(mit);
}

void L1LossMapFun::move_right_if_zero(Coefs &mit){
  move_if_zero(&L1LossMapFun::move_right, mit);
}
void L1LossMapFun::move_left_if_zero(Coefs &mit){
  move_if_zero(&L1LossMapFun::move_left, mit);
}
void L1LossMapFun::move_if_zero(move_it_fun_ptr move, Coefs &mit){
  if(mit.it != loss_map.end() && mit.it->second==0){
    (this->*move)(mit);
    if(mit.it->second==0){
      // if after moving we are still at a zero, they all must be
      // zeros, so move to the end.
      mit.it = loss_map.end();
    }
  }
}

void L1LossMapFun::move_left(Coefs &mit){
  double param_before=get_param(mit);
  double cost_before = get_cost_at_coefs(mit);
  double ldiff = get_Linear_diff(mit);
  if(mit.it == loss_map.begin()){
    mit.it = loss_map.end();
    param_before += MAX_ANGLE;
  }
  mit.it--;
  mit.Linear -= ldiff;
  mit.Constant = cost_before - mit.Linear*param_before;
  moves++;
}

void L1LossMapFun::move_right(Coefs &mit){
  mit.it++;
  if(mit.it == loss_map.end()){
    mit.it = loss_map.begin();
    mit.Constant += mit.Linear*MAX_ANGLE;
  }
  double param_after=get_param(mit);
  double cost_after = mit.Linear*param_after+mit.Constant;
  mit.Linear += get_Linear_diff(mit);
  mit.Constant = cost_after-mit.Linear*param_after;
  moves++;
}

double L1LossMapFun::get_param(L1LossMap::iterator it){
  if(it == loss_map.end())return INFINITY;
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

double L1LossMapFun::get_Linear_diff(Coefs &mit){
  return get_Linear_diff(mit.it);
}

double L1LossMapFun::get_Linear_diff(L1LossMap::iterator it){
  if(it == loss_map.end())return INFINITY;
  return it->second;
}

void L1LossMapFun::push_constant(Cluster new_cl){
  // iterators should be set prior.
  new_cl.sign  = -1;
  new_cl.data_i = data_i;
  new_cl.opt.Linear = new_cl.first.Linear = new_cl.last.Linear = 0;
  new_cl.opt.Constant = new_cl.first.Constant = new_cl.last.Constant = Constant;
  push_cluster(new_cl);
}

void L1LossMapFun::min_with_constant(double constant){
  new_list.clear();
  Constant = constant;
  for
    (ClusterList::iterator it=ptr_list.begin();
     it != ptr_list.end();
     it++){
    double first_cost = get_cost_at_coefs(it->first);
    double opt_cost = get_cost_at_coefs(it->opt);
    double last_cost = get_cost_at_coefs(it->last);
    printf("first=%f opt=%f last=%f constant=%f\n", first_cost, opt_cost, last_cost, constant);
    //first handle this cluster.
    if(first_cost < constant && opt_cost < constant && last_cost < constant){
      printf("case 1 cluster completely below constant push_cluster\n");
      push_cluster(*it);
    }
    if(constant < first_cost && constant < opt_cost && constant < last_cost){
      printf("case 2 cluster completely above constant push_constant\n");
      push_constant(*it);
    }
    if(first_cost < constant && opt_cost < constant && constant < last_cost){
      printf("case 3: convex cluster with a crossing point between opt and last\n");
      Cluster new_cl = *it;
      CrossInfo cinfo = crossing_before(new_cl.last);
      // First push convex piece.
      new_cl.last = cinfo.before;
      push_cluster(new_cl);
      // Then push constant/concave piece.
      Coefs coefs;
      coefs.Constant = constant;
      coefs.Linear = 0;
      new_cl.first = new_cl.opt = new_cl.last = coefs;
      push_cluster(new_cl);
    }
    if(first_cost < constant && constant < opt_cost && constant < last_cost){
      printf("case 4: concave cluster with a crossing point between first and opt\n");
      Cluster new_cl = *it;
      CrossInfo cinfo = crossing_before(new_cl.last);
    }
    if(constant < first_cost && constant < opt_cost && last_cost < constant){
      printf("case 5: concave cluster with a crossing point between opt and last\n");
    }
    if(constant < first_cost && opt_cost < constant && last_cost < constant){
      printf("case 6: cluster with a crossing point between first and opt\n");
    }
    if(first_cost < constant && constant < opt_cost && last_cost < constant){
      printf("case 7: concave cluster with two crossing points\n");
    }
    if(constant < first_cost && opt_cost < constant && constant < last_cost){
      printf("case 8: convex cluster with two crossing points\n");
    }
    // then handle crossing point between this cluster and next.
    ClusterList::iterator next_it=it;
    next_it++;
    bool next_is_first=false;
    if(next_it==ptr_list.end()){
      next_it=ptr_list.begin();
      next_is_first=true;
    }
    double next_first_cost = get_cost_at_coefs(next_it->first);
    if(cost_between(last_cost, constant, next_first_cost)){
      CrossInfo cinfo = crossing_before(next_it->first);
      // TODO handle equality / existing break.
      int new_diff_sign = (last_cost<constant) ? -1 : 1;
      printf("it first=%f opt=%f last=%f next first=%f opt=%f last=%f\n", it->first.Linear, it->opt.Linear, it->last.Linear, next_it->first.Linear, next_it->opt.Linear, it->last.Linear);
      double new_diff = new_diff_sign * it->first.Linear;
      std::pair<L1LossMap::iterator, bool> result;
      result = loss_map.insert(std::pair<double,double>(cinfo.param, new_diff));
      printf("param=%f new_diff=%f result.second=%d\n", cinfo.param, new_diff, result.second);
      L1LossMap::iterator insert_it = result.first;
      Cluster new_cl = *it;
      //cluster before new break.
      if(!(last_cost<constant && it->sign==1)){
	new_cl.last.it = insert_it;
      }
      if(last_cost<constant){
	printf("push_cluster before\n");
	push_cluster(new_cl);
      }else{
	printf("push_constant before\n");
	push_constant(new_cl);
      }
      new_cl = *next_it;
      //cluster after new break.
      if(!(constant < last_cost && it->sign == -1)){
	new_cl.first.it = insert_it;
      }
      if(last_cost<constant){
	printf("push_constant after\n");
	push_constant(new_cl);
      }else{
	printf("push_cluster after\n");
	push_cluster(new_cl);
      }
      if(next_is_first){
	new_list.erase(new_list.begin());
      }
    }
  }
  // loop over new clusters, erase breakpoints in new constant clusters.
  for
    (ClusterList::iterator it=new_list.begin();
     it != new_list.end();
     it++){
    //printf("data_i=%d it->data_i=%d\n", data_i, it->data_i);
    if(it->data_i==data_i){
      //printf("deleting\n");
      Coefs before, after;
      it->opt = it->first;
      after = it->first;
      move_right(after);
      while(after.it != it->last.it){
	before = after;
	move_right(after);
	loss_map.erase(before.it);
      }
    }
  }
  ptr_list = new_list;
}

CrossInfo L1LossMapFun::crossing_before(Coefs coefs){
  double orig_diff = get_cost_at_coefs(coefs)-Constant;
  int orig_sign = sgn(orig_diff);
  int new_sign = orig_sign;
  CrossInfo cinfo;
  // this code works for both kinds of crossing (increasing and decreasing).x
  while(orig_sign == new_sign){
    cinfo.after = coefs;
    move_left(coefs);
    double new_diff = get_cost_at_coefs(coefs)-Constant;
    new_sign = sgn(new_diff);
  }
  cinfo.before = coefs;
  cinfo.param = (Constant-coefs.Constant)/coefs.Linear;
  if(cinfo.param>=MAX_ANGLE){
    cinfo.param -= MAX_ANGLE;
  }
  return cinfo;
}

void L1LossMapFun::push_cluster(const Cluster cl){
  if(new_list.empty()){
    new_list.push_back(cl);
    printf("push first cluster %f %f %f\n", cl.first.it->first, cl.opt.it->first, cl.last.it->first);
    return;
  }
  ClusterList::iterator last_it = new_list.end();
  last_it--;
  if(last_it->sign == cl.sign){
    printf("grow cluster %f -> %f\n", last_it->last.it->first, cl.last.it->first);
    last_it->last = cl.last;
    return;
  }
  if(last_it->last.Linear==0){
    last_it->last.Linear = cl.first.Linear;
  }
  new_list.push_back(cl);
  printf("push new cluster %f %f %f\n", cl.first.it->first, cl.opt.it->first, cl.last.it->first);
}

L1LossMapFun::L1LossMapFun(){
}

Cluster::Cluster(){
}

void Cluster::init(L1LossMap::iterator it_, double Constant_){
}

