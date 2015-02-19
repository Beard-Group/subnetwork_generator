# ifndef IO_H
# define IO_H
#include "../common.h"
using namespace std;
void read_data ( const int N_time_points, const int N_gene, vector <double> &t_d, vector <double> &x_d);
void output_data ( const int gene_ind, const double r0_, const double ea_, const double d_, const string out_r0, const string out_ea, const string out_d, const int start_ka, const int end_ka, const int start_kd, const int end_kd, const vector <int> kavec, const vector <int> kdvec, const host_type &kaval, const host_type &kdval, const string out_nka, const string out_nkd, const string kavec_out_path, const string kdvec_out_path, const string kaval_out_path, const string kdval_out_path);
# endif
