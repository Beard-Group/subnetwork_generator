#ifndef SPLINE_H
#define SPLINE_H
#include "../common.h"
using namespace std;
double pchst ( double arg1, double arg2 );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void spline_pchip_set ( const int n, const vector<double> &x, const vector<double> &f, vector<double> &d );
void spline_deriv ( int gene_index, int N_time_s, vector <double> &t_d_s, vector <double> &x_d_s, vector <double> &deriv_s);
void spline_coeff ( const int N_gene, const int N_time_s, const vector <double> &t_d, const vector <double> &x_d, const vector <double> &deriv_s, vector <double> &cub_coeff_spline);
void calc_spline( const state_type &x_d, const state_type &t_d, const state_type &cub_coeff, const int &ind, state_type &x_spline);
void calc_stats( const state_type &x_spline, state_type &mean, state_type &sd );
struct corr_functor
{
   const double mean_xa, mean_xb;
   corr_functor ( double xa_mean, double xb_mean ) : mean_xa (xa_mean), mean_xb (xb_mean) {}
 
   double operator()(const double &xa_spline, const double &xb_spline)
                    { return ((xa_spline-mean_xa)*(xb_spline-mean_xb));}
};
void calc_corr( const int &N_gene, const state_type &x_spline, const state_type &mean_xd, const state_type &sd_xd, state_type &correlation );
#endif
