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
void calc_stats( const state_type &x_spline, state_type &mean, state_type &sd );
#endif
