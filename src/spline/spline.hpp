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
void calc_spline( const vector<double> &x_d, const vector<double> &t_d, const vector<double> &cub_coeff, const int &ind, vector<double> &x_spline);
void calc_stats( const vector<double> &x_spline, vector<double> &mean, vector<double> &sd );
struct corr_functor
{
   const double mean_xa, mean_xb;
   corr_functor ( double xa_mean, double xb_mean ) : mean_xa (xa_mean), mean_xb (xb_mean) {}
 
   double operator()(const double &xa_spline, const double &xb_spline)
                    { return ((xa_spline-mean_xa)*(xb_spline-mean_xb));}
};
void calc_corr( const int &N_gene, const vector<double> &x_spline, const vector<double> &mean_xd, const vector<double> &sd_xd, vector<double> &correlation );
#endif
