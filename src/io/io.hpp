# ifndef IO_H
# define IO_H
#include "../common.h"
using namespace std;
void read_data ( const int N_time_points, const int N_gene, vector <double> &t_d, vector <double> &x_d);

void output_data ( const int gene_ind, const double r0_, const double ea_, const double d_, const int start_ka, const int end_ka, const int start_kd, const int end_kd, const vector <int> kavec, const vector <int> kdvec, const host_type &kaval, const host_type &kdval, const vector<string>& out_paths);
# endif
