# ifndef SEARCH_H
# define SEARCH_H
# include "../common.h"
# include "../lsoda/cuLsoda.hpp"
# include "../opt/opt.hpp"
/*HOST DEVICE 
unsigned int generate_seed(unsigned int a)
{
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
};*/
struct seed_gen
{
   public:
   //HOST DEVICE 
   unsigned int generate_seed(unsigned int a)
   {
      a = (a+0x7ed55d16) + (a<<12);
      a = (a^0xc761c23c) ^ (a>>19);
      a = (a+0x165667b1) + (a<<5);
      a = (a+0xd3a2646c) ^ (a<<9);
      a = (a+0xfd7046c5) + (a<<3);
      a = (a^0xb55a4f09) ^ (a>>16);
      return a;
   };
   /*HOST DEVICE
   unsigned int operator(const unsigned int ind)()
   {
       return generate_seed(ind);
   }*/
};
void random_conn( const int &gene_ind, unsigned int seed_int, vector <double> &mean, int &n_ka_temp, int &n_kd_temp, vector <int> &kavec_temp, vector <int> &kdvec_temp, vector <double> &kaval_temp, vector <double> &kdval_temp);
/*struct invert_functor
{
   HOST DEVICE
   value_type operator()(const value_type &weight)
   {
       return (1.0 - weight);
   }
};*/
//void biased_select( const int &gene_ind, unsigned int seed_int, vector <double> &mean, const state_type &corr, int &n_ka_temp, int &n_kd_temp, vector <int> &kavec_temp, vector <int> &kdvec_temp, state_type &kaval_temp, state_type &kdval_temp);
# endif 
