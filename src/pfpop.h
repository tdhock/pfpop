#define ERROR_PENALTY_NOT_FINITE 1
#define ERROR_PENALTY_NEGATIVE 2
#define ERROR_DATA_NOT_FINITE 20
#define ERROR_DATA_NEGATIVE 21
#define ERROR_DATA_NOT_LESS_THAN_360 22
#define ERROR_WEIGHT_NOT_FINITE 30
#define ERROR_WEIGHT_NOT_POSITIVE 31

#include <list> //TODO MAP

int pfpop
(const double*, const double, const double*, const int, 
 int*, double*, double*, int*);
int decode
(const int *best_change_ptr,
 const double *best_cost_ptr,
 const double *best_param_ptr,
 const int N_data,
 int *seg_start_ptr,
 int *seg_end_ptr,
 double *seg_param_ptr,
 const int N_segs);

class LinearPiece {
 public:
  double Linear;
  double Constant;
  double min_angle_param;
  double max_angle_param;
  int data_i;
  double prev_angle_param;
  LinearPiece();
  double Loss(double);
  LinearPiece(double li, double co, double m, double M);
  LinearPiece
    (double li, double co, double m, double M, int i, double);
  void print();
};

typedef std::list<LinearPiece> L1LossPieceList;

class PiecewiseLinearLossFun;

typedef void (PiecewiseLinearLossFun::*push_fun_ptr)
(PiecewiseLinearLossFun*,
 PiecewiseLinearLossFun*,
 L1LossPieceList::iterator,
 L1LossPieceList::iterator,
 int);

class PiecewiseLinearLossFun {
 public:
  L1LossPieceList piece_list;
  double weight;
  PiecewiseLinearLossFun();
  void push_sum_pieces(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, L1LossPieceList::iterator, L1LossPieceList::iterator, int);
  void push_min_pieces(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, L1LossPieceList::iterator, L1LossPieceList::iterator, int);
  void while_piece_pairs(PiecewiseLinearLossFun*, PiecewiseLinearLossFun*, push_fun_ptr, int);
  void emplace_piece(double,double,double,double);
  void emplace_piece(double,double,double,double,int,double);
  void enlarge_last_or_emplace(L1LossPieceList::iterator,double,double);
  void enlarge_last_or_emplace(double,double,double,double);
  void enlarge_last_or_emplace(double,double,double,double,int,double);
  void init(double, double);
  void set_to_min_of_one(PiecewiseLinearLossFun *, int);
  void set_to_min_of_two(PiecewiseLinearLossFun *, PiecewiseLinearLossFun *, int);
  void set_to_sum_of(PiecewiseLinearLossFun *, PiecewiseLinearLossFun *, int);
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


