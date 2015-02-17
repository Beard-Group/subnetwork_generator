# ifndef OPT_H
# define OPT_H
# include "../common.h"
# include "../lsoda/cuLsoda.hpp"
double guassrand(); 
void integrate_lsoda_ode ( const int &gene_ind, const vector <double> &x_d, const vector <double> &t_d, const double &mean_xd, const myFex &fex, const myJex &jex, state_type &error_sim_d );
struct error_functor
{
   const value_type m_x, m_s;
   error_functor ( value_type X_m, value_type S_m ) : m_x(X_m), m_s( S_m ) {}
   
   HOST DEVICE
   value_type operator()( const value_type &error_sim, const value_type &y_val)
                        //{ return  error_sim + ((y_val-m_x)/m_s)*((y_val-m_x)/m_s);}
                        { return  error_sim + ((((y_val-m_x)/m_s)*((y_val-m_x)/m_s))/N_time_points);}
};
//
struct update_functor
{
   const value_type m_t_rand;
   value_type temp_p, temp_error, temp_err;
   update_functor ( const value_type t_rand ) : m_t_rand( t_rand ) {}
   // 
   template < class Tuple > 
   HOST DEVICE
   void operator()( Tuple step );
};
struct move_functor
{
   HOST DEVICE 
   unsigned int hash(unsigned int a)
   {
      a = (a+0x7ed55d16) + (a<<12);
      a = (a^0xc761c23c) ^ (a>>19);
      a = (a+0x165667b1) + (a<<5);
      a = (a+0xd3a2646c) ^ (a<<9);
      a = (a+0xfd7046c5) + (a<<3);
      a = (a^0xb55a4f09) ^ (a>>16);
      return a;
   };
   HOST DEVICE
   value_type operator()(unsigned int thread_id, const value_type &param)
   {
       double g;
       clock_t time_seed = clock();
       unsigned int seed = hash(thread_id*time_seed);
       thrust::default_random_engine rng(seed);
       //thrust::uniform_real_distribution<double> u01(-2,2); 
       thrust::random::experimental::normal_distribution<double> u01(1.0,0.5); 
       double r = u01(rng);
       if (r > 1.25) 
       {
          g = r*0.2; 
       } 
       else if (r < 0.75)
       {
          g = r*0.1;
       }
       else g = 0.15;

       return exp10(fmod(log10(param)+g*r,2.0));
   }
};
struct p_functor
{
   const value_type S_MC;
   p_functor ( value_type MC_S ) : S_MC ( MC_S ) {}
   
   HOST DEVICE
   value_type operator()( const value_type &error_opt, const value_type &error_sim)
                        { return  exp((error_opt - error_sim)/S_MC);}
};
void MC_sim ( const int &gene_ind, const vector <double> &x_d, const vector <double> &t_d, const double &mean_xd, const int *n_ka, const int *n_kd, const int &size_ka, const int &size_kd, myFex &fex, const myJex &jex, state_type &error_opt );
# endif
