# ifndef OPT_H
# define OPT_H
# include "../common.h"
# include "../lsoda/cuLsoda.hpp"
double guassrand(); 
void integrate_lsoda_ode ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, const myFex_single &fex, const myJex_single &jex, const int &gene_index,const bool out_val, double &error);
struct move_functor_new
{
   value_type operator()(const value_type &param)
   {
       std::default_random_engine g_gen;
       g_gen.seed(time(NULL));
       std::normal_distribution<double> dist(1.0,0.25);
       double g = 0.2;
       double step_mc = g*dist(g_gen);
       return (fmod(param+step_mc,2.0));
   }
};
void MC_sim_new ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param);
// powof10_functor
struct powof10 : public unary_function<value_type,value_type>
{
       value_type operator()(value_type &x) {return pow(10,x);} 
};
void integrate_lsoda_nm ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_index, const int nka, const int nkd, double &error, const state_type &param);
// log10_functor
struct logb10 : public unary_function<value_type,value_type>
{
       value_type operator()(value_type &x) {return log10(x);} 
};
# endif
