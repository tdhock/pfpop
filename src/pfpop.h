#define ERROR_PENALTY_NOT_FINITE 1
#define ERROR_PENALTY_NEGATIVE 2
#define ERROR_DATA_NOT_FINITE 20
#define ERROR_DATA_NEGATIVE 21
#define ERROR_DATA_NOT_LESS_THAN_360 22
#define ERROR_WEIGHT_NOT_FINITE 30
#define ERROR_WEIGHT_NOT_POSITIVE 31

#include <map>
#include <list>

int pfpop
(const double*, const double, const double*, const int, 
 int*, double*, double*, int*, int*);
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

class LinearCoefsForMap {
 public:
  double Linear;
  double Constant;
  int data_i;
  LinearCoefsForMap();
};

typedef std::list<LinearCoefsForList> L1LossList;
typedef std::map<double, LinearCoefsForMap> L1LossMap;

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


