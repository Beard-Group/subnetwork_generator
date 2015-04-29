# ifndef SEARCH_H
# define SEARCH_H
# include "../common.h"
# include "../lsoda/cuLsoda.hpp"
# include "../optimization/opt.hpp"
void random_select( const int &gene_ind, unsigned int seed_int, const state_type &mean, int &n_ka_temp, int &n_kd_temp, vector <int> &kavec_temp, vector <int> &kdvec_temp, state_type &kaval_temp, state_type &kdval_temp);
# endif 
