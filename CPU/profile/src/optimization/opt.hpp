# ifndef OPT_H
# define OPT_H
# include "../common.h"
# include "../lsoda/cuLsoda.hpp"
double guassrand(); 
void integrate_lsoda_ode ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, const myFex_single &fex, const myJex_single &jex, const int &gene_index,const bool out_val, double &error);
/*struct move_functor
{
   value_type operator()(const value_type &param)
   {
       double g = 0.1;
       return exp10(fmod(log10(param)+g*guassrand(),2.0));
   }
};*/
struct move_functor_new
{
   value_type operator()(const value_type &param)
   {
       double g = 0.0;
       std::default_random_engine g_gen;
       g_gen.seed(time(NULL));
       std::normal_distribution<double> dist(1.0,0.5);
       double r = dist(g_gen);
       if (r > 1.25)
       {
           g = r*0.2;
       }  
       else if (r < 0.75)
       {
           g = r*0.1;
       }
       else g = 0.15;
       double step_mc = g*r;
       return (fmod(param+step_mc,2.0));
   }
};
//void MC_sim ( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex, myJex_single &jex, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param);
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
//void NM_sim( const state_type &x_d, const state_type &t_d, const state_type &mean_xd, myFex_single &fex_nm, myJex_single &jex_nm, const int &gene_ind, const int nka, const int nkd, double &error_ode, state_type &param);
# endif
