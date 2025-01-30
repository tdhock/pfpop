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

class Pointer {
public:
  double Constant, Linear;
  L1LossMap::iterator it;
  Pointer();
  void init(L1LossMap::iterator, double);
};

class L1LossMapFun {
public:
  L1LossMap loss_map;
  Pointer max_ptr, min_ptr;
  Breakpoint last_break;
  double Linear,Constant,min_param,max_param,weight;
  L1LossMapFun();
  void maybe_add(Pointer*);
  double max();
  double argmax();
  double min();
  double argmin();
  void piece(double,double,double,double);
  int  get_data_i(L1LossMap::iterator);
  void set_data_i(L1LossMap::iterator, int);
  double get_Linear_diff(L1LossMap::iterator);
  void   set_Linear_diff(L1LossMap::iterator, double);
  double get_param(L1LossMap::iterator);
  double get_cost_at_pointer(const Pointer);
  Breakpoint* get_break_ptr(L1LossMap::iterator);
  void delete_breaks(double);
  void min_with_constant(double);
  void add_loss_for_data(double,double);
  int move_both_pointers();
  int move_one_pointer(Pointer&, double);
  void move_left(Pointer&);
  void move_right(Pointer&);
  double next_Linear(const Pointer);
  void update_coefs(Pointer&, double, double);
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
  double weight;
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


