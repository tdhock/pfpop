#define pfpop_map_ERROR_PENALTY_NOT_FINITE 1
#define pfpop_map_ERROR_PENALTY_NEGATIVE 2
#define pfpop_map_ERROR_DATA_NOT_FINITE 20
#define pfpop_map_ERROR_DATA_NEGATIVE 21
#define pfpop_map_ERROR_DATA_NOT_LESS_THAN_360 22
#define pfpop_map_ERROR_WEIGHT_NOT_FINITE 30
#define pfpop_map_ERROR_WEIGHT_NOT_POSITIVE 31
#define pfpop_map_ERROR_NOT_IMPLEMENTED 99

#include <map>
#include <list>

int pfpop_map
(const double*, const double, const double*, const int, const char*, 
 double*, double*,
 double*, double*, int*,//max
 double*, double*,
 double*, double*, int*,//min
 int*, int*, int*);

typedef std::map<double, double> L1LossMap;

class Coefs {
public:
  L1LossMap::iterator it;
  double Constant, Linear;
};

class Cluster {
public:
  int sign, data_i;
  Coefs first, last, opt;
  Cluster();
  void init(L1LossMap::iterator, double);
};

class L1LossMapFun;
typedef void (L1LossMapFun::*move_it_fun_ptr) (Coefs&);
typedef std::list<Cluster> ClusterList;

class CrossInfo {
public:
  Coefs before, after;
  double param;
};

class L1LossMapFun {
public:
  L1LossMap loss_map;
  ClusterList ptr_list, new_list;
  ClusterList::iterator cluster_it;
  double Linear,Constant,min_param,max_param,weight,angle;
  double cost;
  int step, moves, data_i;
  L1LossMapFun();
  void all_pointers();
  void update_coefs(Coefs&);
  void move_left(Coefs&);
  void move_left_if_zero(Coefs&);
  void move_right(Coefs&);
  void move_right_if_zero(Coefs&);
  void move_if_zero(move_it_fun_ptr,Coefs&);
  void write_min_or_max(int,int,double*,double*,double*,double*,int*);
  void maybe_move_right(Cluster&,L1LossMap::iterator);
  void maybe_move_erase(L1LossMap::iterator);
  void piece(double,double,double,double);
  void   add_Linear_diff(L1LossMap::iterator, double);
  double prev_Linear(Coefs&);
  double get_Linear_diff(L1LossMap::iterator);
  double get_Linear_diff(Coefs&);
  double get_param(L1LossMap::iterator);
  double get_param(Coefs&);
  void move_to_opt(ClusterList::iterator &it);
  void move_to_diff(L1LossMap::iterator &it, Cluster *p, move_it_fun_ptr);
  double min();
  double max();
  double min_or_max(int);
  void pieces();
  double get_cost_at_coefs(const Coefs);
  void push_cluster(const Cluster);
  void end_move(Cluster&, double);
  double get_param_or_mid(const Cluster);
  void min_with_constant(double);
  void add_loss_for_data(double,double);
  void move_pointers();
  void move_left(L1LossMap::iterator&,Cluster*);
  void move_right(L1LossMap::iterator&,Cluster*);
  CrossInfo crossing_before(Coefs coefs, double constant);
};



