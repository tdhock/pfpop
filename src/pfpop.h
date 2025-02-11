#define ERROR_PENALTY_NOT_FINITE 1
#define ERROR_PENALTY_NEGATIVE 2
#define ERROR_DATA_NOT_FINITE 20
#define ERROR_DATA_NEGATIVE 21
#define ERROR_DATA_NOT_LESS_THAN_360 22
#define ERROR_WEIGHT_NOT_FINITE 30
#define ERROR_WEIGHT_NOT_POSITIVE 31
#define ERROR_NOT_IMPLEMENTED 99

#include <map>
#include <list>

int pfpop_list
(const double*, const double, const double*, const int, 
 int*, double*, double*, int*, int*);
int pfpop_map
(const double*, const double, const double*, const int, 
 int*,
 double*, double*,
 double*, double*,
 double*, double*,
 double*, double*,
 int*, int*, int*);
int decode
(const int *best_change_ptr,
 const double *best_cost_ptr,
 const double *best_param_ptr,
 const int N_data,
 int *seg_start_ptr,
 int *seg_end_ptr,
 double *seg_param_ptr,
 const int N_segs);

class LinearCoefsForList {
 public:
  double Linear;
  double Constant;
  double min_angle_param;
  double max_angle_param;
  double prev_angle_param;
  int data_i;
  LinearCoefsForList();
  double Loss(double);
  LinearCoefsForList(double li, double co, double m, double M);
  LinearCoefsForList
    (double li, double co, double m, double M, int i, double);
  void print();
};

class Breakpoint {
 public:
  double Linear_diff;
  int data_i;
  Breakpoint(double,int);
  Breakpoint();
};

typedef std::list<LinearCoefsForList> L1LossList;
typedef std::map<double, Breakpoint> L1LossMap;


class Mapit {
public:
  L1LossMap::iterator it;
  double get_param();
  void update_coefs(int, double, double,double);
};

class Coefs : public Mapit {
public:
  double Constant, Linear;
  void update_coefs(int,double,double,double);
};

class Cluster {
public:
  Mapit first, last;
  Coefs opt;
  Cluster();
  void init(L1LossMap::iterator, double);
};

class L1LossMapFun;
typedef void (L1LossMapFun::*move_it_fun_ptr) (Mapit&);
typedef std::list<Cluster> ClusterList;

class L1LossMapFun {
public:
  L1LossMap loss_map;
  ClusterList ptr_list;
  Cluster *ptr;
  double Linear,Constant,min_param,max_param,weight,angle;
  int step;
  L1LossMapFun();
  void all_pointers();
  void update_coefs(Mapit&,int,double,double);
  void move_left(Mapit&);
  void move_left_if_zero(Mapit&);
  void move_right(Mapit&);
  void move_right_if_zero(Mapit&);
  void move_if_zero(move_it_fun_ptr,Mapit&);
  void write_min_or_max(int,int,double*,double*,double*,double*);
  void maybe_move_right(Cluster&,L1LossMap::iterator);
  void maybe_move_erase(L1LossMap::iterator);
  void piece(double,double,double,double);
  int  get_data_i(L1LossMap::iterator);
  void set_data_i(L1LossMap::iterator, int);
  double get_Linear_diff(L1LossMap::iterator);
  void   add_Linear_diff(L1LossMap::iterator, double);
  double get_param(L1LossMap::iterator);
  void move_to_diff(L1LossMap::iterator &it, Cluster *p, move_it_fun_ptr);
  double min();
  double max();
  void pieces();
  Cluster* get_min_ptr();
  Cluster* get_max_ptr();
  double get_cost_at_ptr(const Cluster);
  void end_move(Cluster&, double);
  double get_param_or_mid(const Cluster);
  Breakpoint* get_break_ptr(L1LossMap::iterator);
  void delete_breaks(double);
  void min_with_constant(double);
  void add_loss_for_data(double,double);
  void move_pointers();
  void move_left(L1LossMap::iterator&,Cluster*);
  void move_right(L1LossMap::iterator&,Cluster*);
};

class L1LossListFun;

typedef void (L1LossListFun::*push_fun_ptr)
(L1LossListFun*,
 L1LossListFun*,
 L1LossList::iterator,
 L1LossList::iterator,
 int);

class L1LossListFun {
 public:
  L1LossList piece_list;
  double angle, weight;
  int step;
  L1LossListFun();
  void push_sum_pieces(L1LossListFun*, L1LossListFun*, L1LossList::iterator, L1LossList::iterator, int);
  void push_min_pieces(L1LossListFun*, L1LossListFun*, L1LossList::iterator, L1LossList::iterator, int);
  void while_piece_pairs(L1LossListFun*, L1LossListFun*, push_fun_ptr, int);
  void emplace_piece(double,double,double,double);
  void emplace_piece(double,double,double,double,int,double);
  void enlarge_last_or_emplace(L1LossList::iterator,double,double);
  void enlarge_last_or_emplace(double,double,double,double);
  void enlarge_last_or_emplace(double,double,double,double,int,double);
  void init(double, double);
  void set_to_min_of_one(L1LossListFun *, int);
  void set_to_min_of_two(L1LossListFun *, L1LossListFun *, int);
  void set_to_sum_of(L1LossListFun *, L1LossListFun *, int);
  void add(double Constant);
  void multiply(double);
  void print();
  void set_prev_seg_end(int prev_seg_end);
  void findMean(double mean, int *seg_end, double *prev_angle_param);
  double findCost(double mean);
  void Minimize
    (double *best_cost,
     double *best_mean,
     int *data_i);
};


